#include "cooperative_groups.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <array>
#include <assert.h>
#include <atomic>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <stdint.h>
#include <vector>

#include "atomic.cuh"
#include "prefix.cuh"
#include "util.h"
#include "treebuffer.cuh"

template<typename FL, unsigned int levels, unsigned int dimensions> struct GramSchmidtCoeffs
{
private:
  // dimensions x dimensions matrix in row-major storage
  const FL *mu;

public:
  __device__ __host__ inline GramSchmidtCoeffs(const FL *mu)
      : mu(mu)
  {
  }

  __device__ __host__ inline FL get(unsigned int row, unsigned int col)
  {
    assert(row < dimensions);
    assert(col < dimensions);
    return mu[row * dimensions + col];
  }
};

/**
 * Calculates the next child in the enumeration tree for the given node. If the projection of this
 * child is still within the radius, it will be added to the next level in the enumeration tree,
 * with all values correctly initialized.
 *
 * The node is given as a point in the buffer, and if the projection of its child is not within the
 * radius, the subtree of the node is seen as finished and the node is deleted from the buffer.
 *
 * As in all cuda functions, mu should be the triangular basis matrix in column-major storage.
 */
template <typename CG, unsigned int block_size, typename FL, unsigned int levels,
          unsigned int max_nodes_per_level, unsigned int dimensions>
__device__ __host__ void
enumerate_children(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter,
                   OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level, dimensions> &buffer,
                   unsigned int level, GramSchmidtCoeffs<FL, levels, dimensions>& mu, const FL radius_squared,
                   unsigned long long *visited_node_counter)
{

  const unsigned int process_node_count =
      min(buffer.get_node_count(level), cooperative_group.size());
  const bool active = cooperative_group.thread_rank() < process_node_count;
  const unsigned int node_index =
      buffer.get_node_count(level) - process_node_count + cooperative_group.thread_rank();

  if (cooperative_group.thread_rank() == 0)
  {
    atomic_add(visited_node_counter, process_node_count);
  }

  const unsigned int basis_vector_index = levels - level - 1 - 1;

  bool is_children_in_bounds = false;
  FL projected_child_norm_square;
  int child_coefficient = 0;

  unsigned long long start = time();
  if (active)
  {
    CenterToOutIterator &children_iter = buffer.get_children_iterator(level, node_index);
    int next_coefficient_centered      = *children_iter;
    ++children_iter;

    const FL basis_vector_orthogonal_component = mu.get(basis_vector_index, basis_vector_index);

    const FL center_length = buffer.get_center_partsum(level, node_index, basis_vector_index);
    const FL center        = center_length / basis_vector_orthogonal_component;

    child_coefficient = next_coefficient_centered + round(center);

    // the length of the orthogonal part added to the projected point; Equal to alpha
    // times the length of the orthogonalized basis vector
    FL alpha_length = child_coefficient * basis_vector_orthogonal_component - center_length;

    projected_child_norm_square =
        alpha_length * alpha_length + buffer.get_partdist(level, node_index);
    is_children_in_bounds = (projected_child_norm_square <= radius_squared);
  }
  unsigned long long end = time();

  if (cooperative_group.thread_rank() == 0 && active)
  {
    atomic_add(&visited_node_counter[1], end - start);
  }

  start = time();

  const unsigned int child_index =
      buffer.add_open_task(cooperative_group, prefix_counter, is_children_in_bounds, level + 1,
                           child_coefficient, node_index);

  end = time();
  if (cooperative_group.thread_rank() == 0 && active)
  {
    atomic_add(&visited_node_counter[2], end - start);
  }

  start = time();
  if (is_children_in_bounds)
  {
    buffer.set_partdist(level + 1, child_index, projected_child_norm_square);
    // we make a different approach for calculating the center_partsums than other fplll algorithms;
    // instead of calculating the partsums that are necessary for the next step, calculate all
    // partsums that can be calculated with the given information
    for (unsigned int j = 0; j < basis_vector_index; ++j)
    {
      const FL current_center_partsum = buffer.get_center_partsum(level, node_index, j);
      buffer.set_center_partsum(level + 1, child_index, j,
                                current_center_partsum -
                                    child_coefficient * mu.get(j, basis_vector_index));
    }
  }
  end = time();
  if (cooperative_group.thread_rank() == 0 && active)
  {
    atomic_add(&visited_node_counter[3], end - start);
  }

  start = time();
  buffer.filter_open_nodes_segment(cooperative_group, prefix_counter, level, node_index,
                                   is_children_in_bounds, process_node_count);
  end = time();
  if (cooperative_group.thread_rank() == 0 && active)
  {
    atomic_add(&visited_node_counter[4], end - start);
  }
}

template <typename CG, unsigned int block_size, typename FL, unsigned int levels,
          unsigned int max_nodes_per_level, unsigned int dimensions>
__device__ __host__ void init_starting_points(
    CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter, unsigned int group_id,
    unsigned int points_per_group, unsigned int max_start_point_count,
    OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level, dimensions> &buffer, const FL *mu)
{

  const int this_group_start_point_index_begin = group_id * points_per_group;
  const int coefficient = cooperative_group.thread_rank() + this_group_start_point_index_begin;
  const bool active =
      coefficient < max_start_point_count && cooperative_group.thread_rank() < points_per_group;
  const int node_index =
      buffer.add_open_task(cooperative_group, prefix_counter, active, 0, coefficient, 0);

  if (active)
  {
    buffer.set_partdist(0, node_index,
                                           coefficient * coefficient * mu[0] * mu[0]);
    for (unsigned int j = 0; j < levels; ++j)
    {
      buffer.set_center_partsum(0, node_index, j, coefficient * mu[(levels - 1) * levels + j]);
    }
  }
}

/**
 * Processes a found point in the lattice, potentially updating the current radius bound
 * and storing it in the array of possible solutions.
 * 
 * The lattice point is given as a node on the last level of the buffer. Processed points
 * will be deleted from the buffer.
 * 
 * Parameters:
 * shortest_points_per_thread   2d array of shape [total_thread_count, dimension] containing
 *                              the best solution found by each thread until now.
 */
template <typename CG, unsigned int block_size, typename FL, unsigned int levels,
          unsigned int max_nodes_per_level, unsigned int dimension>
__device__ __host__ void
process_solution(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter,
                 OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level, dimension> &buffer,
                 uint32_t *radius_squared_location, FL *shortest_points_per_thread)
{

  const unsigned int process_node_count =
      min(buffer.get_node_count(levels - 1), cooperative_group.size());
  const unsigned int node_index =
      buffer.get_node_count(levels - 1) - process_node_count + cooperative_group.thread_rank();
  bool active = node_index < process_node_count;

  const unsigned int thread_index = thread_id();

  if (active)
  {
    const FL length_squared = buffer.get_partdist(levels - 1, node_index);

    if (abs(length_squared) > .5)
    {
      bool is_shortest_found = length_squared < atomic_min(radius_squared_location, length_squared);
      if (is_shortest_found)
      {
        for (unsigned int i = 0; i < dimension; ++i)
        {
          const FL coefficient = buffer.get_coefficient(levels - 1, node_index, i);
          shortest_points_per_thread[levels * thread_index + i] = coefficient;
        }
      }
    }
  }

  buffer.filter_open_nodes_segment(cooperative_group, prefix_counter, levels - 1, node_index, false,
                                   process_node_count);
}

/**
 * Updates the current tree level on which to work, according to the following three criteria:
 * - The level on which to work should have at least cooperative_group.size() nodes (or as many as 
 * possible), to allow for maximal parallelization. 
 * - Additionally, new nodes will be generated in the next level, and the node buffer for the next 
 * level should have enough capacity.
 * - Last, the lower tree levels should be processed as early as possible, to find solutions quickly,
 * as this can decrease the radius bound and therefore cut tree branches early.
 * 
 * Parameters:
 * current_level        Level on which work was done, will be updated to the level on which to work
 *                      in the next step
 * 
 * cleared_level_count  Count of consecutive levels from the root which have no nodes left (and will
 *                      never get new nodes, as nodes are only generated as children for existing nodes).
 *                      If this is equal to levels, the search is complete
 */
template <unsigned int block_size, typename FL, unsigned int levels,
          unsigned int max_nodes_per_level, unsigned int dimension>
__device__ __host__ void update_current_level(
    OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level, dimension> &buffer,
                     unsigned int group_size, unsigned int &current_level,
                     unsigned int &cleared_level_count)
{
  assert(current_level >= cleared_level_count);

  // if all levels above the processed level are cleared, and in this level all points have been
  // processed, update cleared_levels there is the case that all points have been processed, but no
  // new points have been generated -> possible to clear more than one level
  if (cleared_level_count == current_level)
  {
    while (cleared_level_count < levels && buffer.get_node_count(cleared_level_count) == 0)
    {
      ++cleared_level_count;
    }
  }
  if (cleared_level_count == levels)
  {
    current_level = levels;
    return;
  }

  if (current_level + 1 < levels)
  {
    current_level = max(current_level + 1, cleared_level_count);
  }
  else
  {
    current_level = max(current_level, cleared_level_count);
  }

  // now go up until we find a level with group_size tasks to do; if this does not exist, we stop at
  // the topmost level that is not cleared
  while (current_level > cleared_level_count && buffer.get_node_count(current_level) < group_size)
  {
    --current_level;
  }
  assert(current_level >= cleared_level_count);
  assert(current_level == levels || buffer.get_node_count(current_level) > 0);
}

/**
 * Searches the expansion tree for a shortest nonzero lattice vector. This vector will
 * be written into any entry of shortest_points_per_thread.
 * 
 * Parameters:
 * radius_squared_location      Location accessible by all threads that will be atomically updated with
 *                              the squared norm of the shortest found (nonzero) lattice vector
 * 
 * shortest_points_per_thread   2d array of shape [total_thread_count, levels] that will be filled
 *                              with short found vectors
 */
template <typename CG, unsigned int block_size, typename FL,
          unsigned int levels, unsigned int max_nodes_per_level, unsigned int dimensions>
__device__ __host__ void
search_enumeration_tree(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter,
                   OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level, dimensions> &buffer,
                   GramSchmidtCoeffs<FL, levels, dimensions>& mu, FL *shortest_points_per_thread, uint32_t* radius_squared_location,
                   unsigned long long *visited_node_counter)
{

  unsigned int level               = 0;
  unsigned int cleared_level_count = 0;

  while (cleared_level_count < levels)
  {
    assert(buffer.get_node_count(level) > 0);
    assert(level < levels);

    if (level + 1 < levels)
    {
      enumerate_children<CG, block_size, FL, levels, max_nodes_per_level, dimensions>(
          cooperative_group, prefix_counter, buffer, level, mu, int_to_float_order_preserving_bijection(*radius_squared_location),
          visited_node_counter);
    }
    else
    {
      process_solution<CG, block_size, FL, levels, max_nodes_per_level, dimensions>(
          cooperative_group, prefix_counter, buffer, radius_squared_location,
          shortest_points_per_thread);
    }

    if (TRACE)
    {
      cooperative_group.sync();
      debug_message_thread("Worked on level %d\n", level);
      for (unsigned int i = (unsigned int)max((int)level - 2, 0); i < min(level + 3, levels); ++i)
      {
        debug_message_thread("  %d: %d\n", i, buffer.get_node_count(i));
      }
    }

    // sync to allow access to the count of nodes on the current level
    cooperative_group.sync();
    update_current_level<block_size, FL, levels, max_nodes_per_level>(
        buffer, cooperative_group.size(), level, cleared_level_count);
    // sync to allow writes to the count of nodes on the current level
    cooperative_group.sync();
    assert(level <= levels);
  }
}

/**
 * Kernel that wraps the parallel execution of the enumeration tree search.
 * 
 * mu               dimensions x dimensions matrix in column-major storage
 * 
 * starting_points  2d array of shape [starting_point_count, dimensions - levels + 1] of coefficients
 */
template <unsigned int block_size, typename FL, unsigned int levels,
          unsigned int max_nodes_per_level, unsigned int dimensions>
__global__ void __launch_bounds__(512, 1)
    search_enumeration_tree_kernel(const FL *starting_points,
                                   const unsigned int starting_point_count,
                                   const unsigned int starting_points_per_group,
                                   unsigned int* processed_starting_points,
                                   unsigned char *buffer_memory, 
                                   uint32_t *radius_squared_location, 
                                   const FL *gram_schmidt_coeffs,
                                   FL *shortest_points_per_thread,
                                   unsigned long long *visited_node_counter)
{
  static_assert(levels >= 1, "Enumeration tree must have at least a root");

  extern __shared__ unsigned char shared_mem[PrefixCounter<cooperative_groups::thread_block,
                                                           block_size>::shared_mem_size_in_bytes() + sizeof(unsigned int)];

  typedef cooperative_groups::thread_block CG;
  CG group = cooperative_groups::this_thread_block();
  PrefixCounter<CG, block_size> prefix_counter(shared_mem + sizeof(unsigned int));
  const unsigned int group_id = blockIdx.x;

  typedef OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level, dimensions> NodeBuffer;
  NodeBuffer buffer(buffer_memory + group_id * NodeBuffer::memory_size_in_bytes_per_group);

  GramSchmidtCoeffs<FL, levels, dimensions> mu(gram_schmidt_coeffs);

  assert(starting_points_per_group <= group.size());
  unsigned int *shared_starting_point_offset = reinterpret_cast<unsigned int*>(shared_mem);
  if (group.thread_rank() == 0)
  {
    *shared_starting_point_offset =
        atomic_add(processed_starting_points, starting_points_per_group);
  }
  group.sync();
  while (*shared_starting_point_offset < starting_point_count)
  {
    buffer.init(group);
    group.sync();

    const unsigned int starting_point_index = *shared_starting_point_offset + group.thread_rank();
    const bool active                       = group.thread_rank() < starting_points_per_group &&
                        starting_point_index < starting_point_count;

    const unsigned int node_index =
        buffer.add_starting_task(group, prefix_counter, active);

    if (active)
    {
      constexpr unsigned int starting_point_last_dim = dimensions - levels + 1;
      FL partdist                                    = 0;
      for (unsigned int i = levels - 1; i < dimensions; ++i)
      {
        const FL coefficient = starting_points[starting_point_index * starting_point_last_dim + i];
        buffer.set_start_coefficient(node_index, i, coefficient);
        const FL orthogonal_component = mu.get(i, i);
        partdist += coefficient * coefficient * orthogonal_component * orthogonal_component;
      }
      buffer.set_partdist(0, node_index, partdist);
      for (unsigned int basis_vector_index = 0; basis_vector_index < levels; ++basis_vector_index)
      {
        FL center_partsum = 0;
        for (unsigned int i = levels - 1; i < dimensions; ++i)
        {
          const FL coefficient =
              starting_points[starting_point_index * starting_point_last_dim + i];
          center_partsum -= coefficient * mu.get(basis_vector_index, i);
        }
        buffer.set_center_partsum(0, node_index, basis_vector_index, center_partsum);
      }
    }

    group.sync();

    search_enumeration_tree(group, prefix_counter, buffer, mu, shortest_points_per_thread,
                            radius_squared_location, visited_node_counter);

    if (group.thread_rank() == 0)
    {
      *shared_starting_point_offset =
          atomic_add(processed_starting_points, starting_points_per_group);
    }
    group.sync();
  }
}

template <typename FL, unsigned int levels, unsigned int dimensions>
CudaPtr<FL> search_enumeration(const std::array<std::array<FL, dimensions>, dimensions>& host_mu, const std::vector<std::array<FL, dimensions - levels + 1>>& starting_points, unsigned int &output_point_count, const FL initial_radius)
{

  constexpr unsigned int block_size          = 512;
  constexpr unsigned int max_nodes_per_level = 2 * block_size;

  typedef OpenExpansionTreeNodeBuffer<FL, levels, block_size * 2, dimensions> NodeBuffer;

  const uint32_t radius_store_init =
      float_to_int_order_preserving_bijection(initial_radius * initial_radius);
  CudaPtr<uint32_t> radius_squared = alloc(uint32_t, 1);
  check(cudaMemcpy(radius_squared.get(), &radius_store_init, sizeof(uint32_t),
                   cudaMemcpyHostToDevice));

  constexpr unsigned int starting_point_dim    = dimensions - levels + 1;
  CudaPtr<FL> device_starting_points = alloc(FL, starting_points.size() * starting_point_dim);
  for (unsigned int i = 0; i < starting_points.size(); ++i)
  {
    check(cudaMemcpy(&device_starting_points.get()[i * starting_point_dim],
                     starting_points[i].data(), starting_point_dim * sizeof(FL),
                     cudaMemcpyHostToDevice));
  }

  const unsigned int starting_points_per_group = 4;
  const unsigned int block_count               = std::min<unsigned int>(64, (starting_points.size() - 1) / starting_points_per_group + 1);

  CudaPtr<unsigned int> processed_starting_points = alloc(unsigned int, 1);
  check(cudaMemset(processed_starting_points.get(), 0, sizeof(unsigned int)));

  CudaPtr<FL> shortest_points_per_thread = alloc(FL, block_count * block_size * levels);

  CudaPtr<FL> mu = alloc(FL, dimensions * dimensions);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    check(cudaMemcpy(&mu.get()[i * dimensions], host_mu[i].data(), dimensions * sizeof(FL),
                     cudaMemcpyHostToDevice));
  }

  std::cout << "Start " << block_count << " blocks with " << block_size
            << " threads, global memory "
            << (block_count * NodeBuffer::memory_size_in_bytes_per_group) << " bytes" << std::endl;
  CudaPtr<unsigned char> buffer_memory =
      alloc(unsigned char, block_count *NodeBuffer::memory_size_in_bytes_per_group);
  CudaPtr<unsigned long long> visited_node_counter = alloc(unsigned long long, 5);

  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

  search_enumeration_tree_kernel<block_size, FL, levels, max_nodes_per_level, dimensions>
      <<<dim3(block_count), dim3(block_size)>>>(
          device_starting_points.get(), starting_points.size(), starting_points_per_group,
          processed_starting_points.get(), buffer_memory.get(), radius_squared.get(), mu.get(),
          shortest_points_per_thread.get(), visited_node_counter.get()
      );
  check(cudaDeviceSynchronize());
  check(cudaGetLastError());

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  uint32_t end_radius_square;
  check(cudaMemcpy(&end_radius_square, radius_squared.get(), sizeof(uint32_t),
                   cudaMemcpyDeviceToHost));
  std::cout << "End radius: " << sqrt(int_to_float_order_preserving_bijection(end_radius_square))
            << std::endl;

  unsigned long long counter;
  check(cudaMemcpy(&counter, visited_node_counter.get(), sizeof(unsigned long long),
                   cudaMemcpyDeviceToHost));
  std::cout << "Visited " << counter << " nodes in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
            << std::endl;

  print_performance_counter(&visited_node_counter.get()[1]);
  print_performance_counter(&visited_node_counter.get()[2]);
  print_performance_counter(&visited_node_counter.get()[3]);
  print_performance_counter(&visited_node_counter.get()[4]);

  output_point_count = block_size * block_count;
  return shortest_points_per_thread;
}

/*
template <typename FL, unsigned int min_levels, unsigned int max_levels>
CudaPtr<FL> search_enumeration_choose_template_instance(const FL *mu, const unsigned int levels,
                                                        const FL initial_radius)
{
  if constexpr (min_levels < max_levels)
  {
    constexpr unsigned int mid = (min_levels + max_levels) / 2;
    if (levels <= mid)
    {
      return search_enumeration_choose_template_instance<FL, min_levels, mid>(mu, levels,
                                                                              initial_radius);
    }
    else
    {
      return search_enumeration_choose_template_instance<FL, mid + 1, max_levels>(mu, levels,
                                                                                  initial_radius);
    }
  }
  else
  {
    static_assert(min_levels == max_levels,
                  "Error in template code, expected min_levels >= max_levels");
    return search_enumeration<FL, min_levels>(mu, levels);
  }
}

constexpr unsigned int max_levels_for_enumeration = 128;

std::unique_ptr<float[]> search_enumeration_cuda(const float *mu, const unsigned int levels, const FL initial_radius)
{
  assert(levels <= max_levels_for_enumeration);
  CudaPtr<float> result =
      search_enumeration_choose_template_instance<float, 1, max_levels_for_enumeration>(mu, levels, initial_radius);
  std::unique_ptr<float[]> host_result(new float[levels]);
  check(
      cudaMemcpy(host_result.get(), result.get(), levels * sizeof(float), cudaMemcpyDeviceToHost));
  return host_result;
}
*/