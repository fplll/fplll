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

#include "atomic.cuh"
#include "prefix.cuh"
#include "testdata.h"
#include "util.h"

constexpr unsigned int EXPAND_NODES_PER_CHILDREN = 1;

/**
Provides an iterator that iterates over Z in the order 0, -1, 1, -2, 2, -3, 3, ...
*/
struct CenterToOutIterator
{
  int current;

  __device__ __host__ inline CenterToOutIterator() : current(0) {}

  __device__ __host__ inline void operator++()
  {
    current += static_cast<int>(current >= 0);
    current = -current;
  }

  __device__ __host__ inline void operator+=(unsigned int val)
  {
    if (val % 2 == 0)
    {
      current += current >= 0 ? val / 2 : -val / 2;
    }
    else
    {
      operator+=(val - 1);
      operator++();
    }
  }

  __device__ __host__ inline int operator*() { return current; }
};

/**
 * Stores a set of enumeration tree nodes whose subtrees have to be processed, ordered by tree level in a way
 * that maximizes memory coalescing.
 * 
 * For general information about how the tree search works, refer to the readme.
 */
template <typename FL, unsigned int levels, unsigned int max_nodes_per_level>
struct OpenExpansionTreeNodeBuffer
{
  // shape [levels, levels, max_nodes_per_level]; origin is not contained (would be at level -1)
  FL *__restrict__ coefficients;
  // shape [levels, max_nodes_per_levle]
  CenterToOutIterator *__restrict__ children_iterator;

  // shape [levels, max_nodes_per_level], contains the squared norm of the point projected into the
  // subspace perpendicular to the first k lattice vectors
  FL *__restrict__ partdist;
  // shape [levels, levels, max_nodes_per_level]; [k, j, i] contains the scalar product of the i-th point
  // at level k with the j-th orthogonalized basis vector
  FL *__restrict__ center_partsums;

  // shape [levels]
  unsigned int *__restrict__ open_node_count;

  __device__ __host__ OpenExpansionTreeNodeBuffer(unsigned char *memory)
      : coefficients(reinterpret_cast<FL *>(memory)),
        children_iterator(reinterpret_cast<CenterToOutIterator *>(
            memory + (levels * sizeof(FL)) * levels * max_nodes_per_level)),
        partdist(
            reinterpret_cast<FL *>(memory + (levels * sizeof(FL) + sizeof(CenterToOutIterator)) *
                                                levels * max_nodes_per_level)),
        center_partsums(reinterpret_cast<FL *>(
            memory + (levels * sizeof(FL) + sizeof(CenterToOutIterator) + sizeof(FL)) * levels *
                         max_nodes_per_level)),
        open_node_count(reinterpret_cast<unsigned int *>(
            memory +
            (levels * sizeof(FL) + sizeof(CenterToOutIterator) + sizeof(FL) + levels * sizeof(FL)) *
                levels * max_nodes_per_level))
  {
  }

  template <typename CG> __device__ __host__ inline void init(CG &cooperative_group)
  {
    if (cooperative_group.thread_rank() == 0)
    {
      for (unsigned int i = 0; i < levels; ++i)
      {
        open_node_count[i] = 0;
      }
    }
  }

  constexpr static size_t memory_size_in_bytes_per_group =
      (levels * sizeof(int) + sizeof(CenterToOutIterator) + sizeof(FL) + levels * sizeof(FL)) *
          levels * max_nodes_per_level +
      sizeof(unsigned int) * levels;

  static_assert(memory_size_in_bytes_per_group < std::numeric_limits<unsigned int>::max(),
                "Requires more memory than indexable with unsigned int");

  __device__ __host__ inline unsigned int get_node_count(unsigned int tree_level)
  {
    return open_node_count[tree_level];
  }

  /**
   * Returns the squared norm of the node_index-th point at level tree_level projected into
   * the subspace perpendicular to the first levels - tree_level basis vectors
   */
  __device__ __host__ inline FL get_partdist(unsigned int tree_level,
                                                                unsigned int node_index)
  {
    return partdist[tree_level * max_nodes_per_level + node_index];
  }

  /**
   * Sets the squared norm of the node_index-th point at level tree_level projected into
   * the subspace perpendicular to the first levels - tree_level basis vectors
   */
  __device__ __host__ inline void set_partdist(unsigned int tree_level,
                                                                  unsigned int node_index, FL value)
  {
    partdist[tree_level * max_nodes_per_level + node_index] = value;
  }

  /**
   * Returns the inner product of the node_index-th point at level tree_level with
   * the given orthogonalized basis vector (given via index in basis)
   */
  __device__ __host__ inline FL get_center_partsum(unsigned int tree_level, unsigned int node_index,
                                                  unsigned int orthogonalized_basis_vector_index)
  {
    return center_partsums[tree_level * levels * max_nodes_per_level +
                       orthogonalized_basis_vector_index * max_nodes_per_level + node_index];
  }

  /**
   * Returns the inner product of the node_index-th point at level tree_level with
   * the given orthogonalized basis vector (given via index in basis)
   */
  __device__ __host__ inline void set_center_partsum(unsigned int tree_level,
                                                    unsigned int node_index,
                                                    unsigned int orthogonalized_basis_vector_index,
                                                    FL value)
  {
    center_partsums[tree_level * levels * max_nodes_per_level +
                orthogonalized_basis_vector_index * max_nodes_per_level + node_index] = value;
  }

  /**
   * Returns the iterator over the children coefficients
   */
  __device__ __host__ inline CenterToOutIterator &get_children_iterator(unsigned int tree_level,
                                                                        unsigned int node_index)
  {
    return children_iterator[tree_level * max_nodes_per_level + node_index];
  }

  /**
   * Returns the coefficient of the node_index-th point at level tree_level w.r.t.
   * the given basis vector (given by index). Called x in other fplll algorithms
   */
  __device__ __host__ inline FL get_coefficient(unsigned int tree_level, unsigned int node_index,
                                                unsigned int coordinate)
  {
    return coefficients[tree_level * max_nodes_per_level * levels +
                        coordinate * max_nodes_per_level + node_index];
  }

  /**
   * Adds a new node to the enumeration tree at the given tree level. The coefficient of this
   * node must be known, but the additionally calculated data (center_partsums, partdist) must
   * be set later via the corresponding setter functions.
   * 
   * This modifies the node count and the nodes at the given tree level, so any subsequent
   * reads or writes on this data must be synchronized.
   * 
   * Returns the index of the created node.
   */
  template <typename CG, unsigned int block_size>
  __device__ __host__ inline unsigned int
  add_open_task(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter,
                bool this_thread_active, unsigned int tree_level, FL coefficient,
                unsigned int parent_node_index)
  {
    unsigned int total_new_tasks      = 0;
    const unsigned int existing_tasks = get_node_count(tree_level);
    const unsigned int new_task_offset =
        prefix_counter.prefix_count(cooperative_group, this_thread_active, total_new_tasks);
    assert(existing_tasks + total_new_tasks <= max_nodes_per_level);
    const unsigned int new_task_index = new_task_offset + existing_tasks;

    if (this_thread_active)
    {
      assert(tree_level == 0 || parent_node_index < open_node_count[tree_level - 1]);
      children_iterator[tree_level * max_nodes_per_level + new_task_index] = CenterToOutIterator{};
      coefficients[tree_level * max_nodes_per_level * levels +
                   (levels - tree_level - 1) * max_nodes_per_level + new_task_index] = coefficient;
      for (unsigned int i = levels - tree_level; i < levels; ++i)
      {
        coefficients[tree_level * max_nodes_per_level * levels + i * max_nodes_per_level +
                     new_task_index] = get_coefficient(tree_level - 1, parent_node_index, i);
      }
    }

    if (cooperative_group.thread_rank() == 0)
    {
      assert(total_new_tasks + open_node_count[tree_level] <= max_nodes_per_level);
      open_node_count[tree_level] += total_new_tasks;
    }
    return new_task_index;
  }

  /**
   * Filters the open nodes at the given tree level. All nodes for which a thread invokes this function
   * with parameter keep_this_thread_task == false are removed, all other tasks on the tree level are
   * kept.

   * To allow an efficient implementation, it is required that node_index ==
   * open_nodes_count[tree_level] - active_thread_count + cooperative_group.thread_rank(); For inactive
   * threads (i.e. cooperative_group.thread_rank() >= active_thread_count), the value for
   * keep_this_thread_task is ignored, as there is no task to remove.
   * 
   * This modifies the node count and the nodes at the given tree level, so any subsequent
   * reads or writes on this data must be synchronized.
   */
  template <typename CG, unsigned int block_size>
  __device__ __host__ inline void
  filter_open_nodes_segment(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter,
                            unsigned int tree_level, unsigned int node_index,
                            bool keep_this_thread_task, unsigned int active_thread_count)
  {
    assert(active_thread_count <= open_node_count[tree_level]);
    assert(node_index == open_node_count[tree_level] - active_thread_count +
                             cooperative_group.thread_rank());

    unsigned int kept_tasks = 0;
    const bool is_active =
        keep_this_thread_task && cooperative_group.thread_rank() < active_thread_count;
    const unsigned int new_offset =
        prefix_counter.prefix_count(cooperative_group, is_active, kept_tasks);
    const unsigned int new_index =
        new_offset + open_node_count[tree_level] - active_thread_count;

    const unsigned int old_node_index = tree_level * max_nodes_per_level + node_index;
    const unsigned int new_node_index = tree_level * max_nodes_per_level + new_index;

    FL projected_point_norm_square_tmp;
    FL partial_sum_tmp[levels];
    int coefficient_tmp[levels];
    CenterToOutIterator children_iterator_tmp;
    if (is_active)
    {
      projected_point_norm_square_tmp = partdist[old_node_index];
      children_iterator_tmp           = children_iterator[old_node_index];
      for (unsigned int j = 0; j < levels; ++j)
      {
        partial_sum_tmp[j] = center_partsums[tree_level * levels * max_nodes_per_level +
                                         j * max_nodes_per_level + node_index];
      }
      for (unsigned int j = levels - tree_level - 1; j < levels; ++j)
      {
        coefficient_tmp[j] = coefficients[tree_level * levels * max_nodes_per_level +
                                          j * max_nodes_per_level + node_index];
      }
    }

    cooperative_group.sync();

    if (is_active)
    {
      partdist[new_node_index] = projected_point_norm_square_tmp;
      children_iterator[new_node_index]            = children_iterator_tmp;
      for (unsigned int j = 0; j < levels; ++j)
      {
        center_partsums[tree_level * levels * max_nodes_per_level + j * max_nodes_per_level +
                    new_index] = partial_sum_tmp[j];
      }
      for (unsigned int j = levels - tree_level - 1; j < levels; ++j)
      {
        coefficients[tree_level * levels * max_nodes_per_level + j * max_nodes_per_level +
                     new_index] = coefficient_tmp[j];
      }
    }

    if (cooperative_group.thread_rank() == 0)
    {
      open_node_count[tree_level] -= active_thread_count - kept_tasks;
    }
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
          unsigned int max_nodes_per_level>
__device__ __host__ void
enumerate_children(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter,
                   OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level> &buffer,
                   unsigned int level, const FL *mu, const FL radius_squared,
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

  if (active)
  {
    CenterToOutIterator &children_iter = buffer.get_children_iterator(level, node_index);
    int next_coefficient_centered      = *children_iter;
    ++children_iter;

    const FL basis_vector_orthogonal_component =
        mu[basis_vector_index * levels + basis_vector_index];

    const FL center_length = buffer.get_center_partsum(level, node_index, basis_vector_index);
    const FL center        = round(center_length / basis_vector_orthogonal_component);

    child_coefficient = next_coefficient_centered + center;

    // the length of the orthogonal part added to the projected point; Equal to alpha
    // times the length of the orthogonalized basis vector
    FL alpha_length = child_coefficient * basis_vector_orthogonal_component - center_length;

    projected_child_norm_square = alpha_length * alpha_length +
                                  buffer.get_partdist(level, node_index);
    is_children_in_bounds = (projected_child_norm_square <= radius_squared);
  }

  const unsigned int child_index =
      buffer.add_open_task(cooperative_group, prefix_counter, is_children_in_bounds, level + 1,
                           child_coefficient, node_index);

  if (is_children_in_bounds)
  {
    buffer.set_partdist(level + 1, child_index, projected_child_norm_square);
    // we make a different approach for calculating the center_partsums than other fplll algorithms; 
    // instead of calculating the partsums that are necessary for the next step, calculate all partsums 
    // that can be calculated with the given information
    for (unsigned int j = 0; j < basis_vector_index; ++j)
    {
      const FL current_center_partsum = buffer.get_center_partsum(level, node_index, j);
      buffer.set_center_partsum(
          level + 1, child_index, j,
          current_center_partsum -
              child_coefficient * mu[basis_vector_index * levels + j]);
    }
  }

  buffer.filter_open_nodes_segment(cooperative_group, prefix_counter, level, node_index,
                                   is_children_in_bounds, process_node_count);
}

/**
 * Calculates the count of multiples of the last vector so that its projection onto the
 * last dimension is less or equal than radius. Just returns the count of nonzero multiples.
 * 
 * As in all cuda functions, mu should be the triangular basis matrix in column-major storage.
 */
template <typename FL, unsigned int levels>
constexpr unsigned int calc_starting_point_count(const FL *mu, const FL radius)
{
  return static_cast<unsigned int>(floor(abs(radius / mu[levels * (levels - 1) + levels - 1]))) + 1;
}

template <typename CG, unsigned int block_size, typename FL, unsigned int levels,
          unsigned int max_nodes_per_level>
__device__ __host__ void init_starting_points(
    CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter, unsigned int group_id,
    unsigned int points_per_group, unsigned int max_start_point_count,
    OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level> &buffer, const FL *mu)
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
 * shortest_points_per_thread   2d array of shape [total_thread_count, levels] containing
 *                              the best solution found by each thread until now.
 */
template <typename CG, unsigned int block_size, typename FL,
          unsigned int levels, unsigned int max_nodes_per_level>
__device__ __host__ void
process_solution(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter,
                 OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level> &buffer,
                 uint32_t *radius_squared_location, FL *shortest_points_per_thread)
{

  const unsigned int process_node_count =
      min(buffer.get_node_count(levels - 1), cooperative_group.size());
  const unsigned int node_index =
      buffer.get_node_count(levels - 1) - process_node_count + cooperative_group.thread_rank();
  bool active = node_index < process_node_count;

  const unsigned int thread_index = thread_id();
  assert(cooperative_group.thread_rank() == thread_index);

  if (active)
  {
    const FL length_squared = buffer.get_partdist(levels - 1, node_index);

    if (abs(length_squared) > .5)
    {
      bool is_shortest_found = length_squared < atomic_min(radius_squared_location, length_squared);
      if (is_shortest_found)
      {
        for (unsigned int i = 0; i < levels; ++i)
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
          unsigned int max_nodes_per_level>
__device__ __host__ void
update_current_level(OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level> &buffer,
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
          unsigned int levels, unsigned int max_nodes_per_level>
__device__ __host__ void
search_enumeration_tree(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter,
                   OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level> &buffer,
                   const FL *mu, FL *shortest_points_per_thread, uint32_t* radius_squared_location,
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
      enumerate_children<CG, block_size, FL, levels, max_nodes_per_level>(
          cooperative_group, prefix_counter, buffer, level, mu, int_to_float_order_preserving_bijection(*radius_squared_location),
          visited_node_counter);
    }
    else
    {
      process_solution<CG, block_size, FL, levels, max_nodes_per_level>(
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
 * TODO: initialization of starting points
 */
template <unsigned int block_size, typename FL, unsigned int levels,
          unsigned int max_nodes_per_level>
__global__ void __launch_bounds__(512, 1)
    search_enumeration_tree_kernel(const unsigned int starting_points_per_group,
                              const unsigned int total_starting_points,
                              unsigned char *buffer_memory, uint32_t *radius_squared_location, const FL *mu,
                              FL *shortest_points_per_thread,
                              unsigned long long *visited_node_counter)
{

  extern __shared__ unsigned char shared_mem[PrefixCounter<cooperative_groups::thread_block,
                                                           block_size>::shared_mem_size_in_bytes()];

  typedef cooperative_groups::thread_block CG;
  CG group = cooperative_groups::this_thread_block();
  PrefixCounter<CG, block_size> prefix_counter(shared_mem);
  const unsigned int group_id = blockIdx.x;

  typedef OpenExpansionTreeNodeBuffer<FL, levels, max_nodes_per_level> NodeBuffer;
  NodeBuffer buffer(buffer_memory + group_id * NodeBuffer::memory_size_in_bytes_per_group);
  buffer.init(group);

  init_starting_points(group, prefix_counter, group_id, starting_points_per_group,
                       total_starting_points, buffer, mu);
  group.sync();
  search_enumeration_tree(group, prefix_counter, buffer, mu, shortest_points_per_thread,
                          radius_squared_location, visited_node_counter);
}

template <typename FL, unsigned int levels>
CudaPtr<FL> search_enumeration(const FL* host_mu, unsigned int &output_point_count, const FL initial_radius)
{

  constexpr unsigned int block_size          = 512;
  constexpr unsigned int max_nodes_per_level = (EXPAND_NODES_PER_CHILDREN + 1) * block_size;

  typedef OpenExpansionTreeNodeBuffer<FL, levels, block_size * 2> NodeBuffer;

  const uint32_t radius_store_init =
      float_to_int_order_preserving_bijection(initial_radius * initial_radius);
  CudaPtr<uint32_t> radius_squared = alloc(uint32_t, 1);
  check(cudaMemcpy(radius_squared.get(), &radius_store_init, sizeof(uint32_t),
                   cudaMemcpyHostToDevice));

  const unsigned int starting_points_per_group = 4;
  const unsigned int total_starting_point_count =
      calc_starting_point_count<FL, levels>(host_mu, initial_radius);
  const unsigned int block_count = (total_starting_point_count - 1) / starting_points_per_group + 1;
  CudaPtr<FL> shortest_points_per_thread = alloc(FL, block_count * block_size * levels);

  CudaPtr<FL> mu = alloc(FL, levels * levels);
  check(cudaMemcpy(mu.get(), host_mu, levels * levels * sizeof(FL), cudaMemcpyHostToDevice));

  std::cout << "Start " << block_count << " blocks with " << block_size
            << " threads, global memory "
            << (block_count * NodeBuffer::memory_size_in_bytes_per_group) << " bytes" << std::endl;
  CudaPtr<unsigned char> buffer_memory =
      alloc(unsigned char, block_count *NodeBuffer::memory_size_in_bytes_per_group);
  CudaPtr<unsigned long long> visited_node_counter = alloc(unsigned long long, 1);

  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

  search_enumeration_tree_kernel<block_size, FL, levels, max_nodes_per_level>
      <<<dim3(block_count), dim3(block_size)>>>(
          starting_points_per_group, total_starting_point_count, buffer_memory.get(),
          radius_squared.get(), mu.get(), shortest_points_per_thread.get(),
          visited_node_counter.get());
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

  output_point_count = block_size * block_count;
  return shortest_points_per_thread;
}

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

// instantiating this greatly increases compile time

/*
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