#ifndef FPLLL_ENUM_CUH
#define FPLLL_ENUM_CUH

#include "cooperative_groups.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <assert.h>
#include <chrono>
#include <iostream>
#include <limits>
#include <stdint.h>
#include <functional>
#include <vector>

#include "atomic.h"
#include "cuda_util.cuh"
#include "prefix.cuh"
#include "streaming.cuh"
#include "recenum.cuh"

namespace cuenum
{

constexpr bool CUENUM_TRACE = false;

/**
 * Stores the state of the enumeration tree search, i.e. all nodes of the tree whose subtrees have
 * to be searched, or are currently searched, ordered by tree level.
 *
 * A node is a dimensions_per_level-level subtree of the enumeration tree that is traversed via
 * enumerate_recursive() and is given by the point corresponding to its root, i.e. a point in the
 * sublattice spanned by the last level * dimensions_per_level basis vectors.
 */
template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct SubtreeEnumerationBuffer
{
  // private:
  // coefficients of the children enumeration for this point, used to pause and resume
  // enumerate_recursive() shape [levels, dimensions_per_level, max_nodes_per_level]
  enumi *enumeration_x;
  // last dimensions_per_level coefficients of the point, the other coefficients must be retrieved
  // by querying the parent node shape [levels, dimensions_per_level, max_nodes_per_level]
  enumi *coefficients;
  // inner products with the scaled lattice vectors of the point, only the first (levels - level) *
  // dimensions_per_level are of interest shape [levels, dimensions, max_nodes_per_level]
  enumf *center_partsum;
  // squared norm of the point, projected into the perpendicular subspace to the first (levels -
  // level) * dimensions_per_level basis vectors shape [levels, max_nodes_per_level]
  enumf *partdist;
  // shape [levels, max_nodes_per_level]
  unsigned int *parent_indices;
  // shape [levels]
  unsigned int *open_node_count;

  constexpr static unsigned int dimensions = levels * dimensions_per_level;

  constexpr static unsigned int enumeration_x_size_in_bytes =
      sizeof(enumi) * levels * dimensions_per_level * max_nodes_per_level;

  constexpr static unsigned int coefficient_size_in_bytes =
      sizeof(enumi) * levels * dimensions_per_level * max_nodes_per_level;

  constexpr static unsigned int center_partsum_size_in_bytes =
      sizeof(enumf) * levels * dimensions * max_nodes_per_level;

  constexpr static unsigned int partdist_size_in_bytes =
      sizeof(enumf) * levels * max_nodes_per_level;

  constexpr static unsigned int parent_indices_size_in_bytes =
      sizeof(unsigned int) * levels * max_nodes_per_level;

  constexpr static unsigned int open_node_count_size_in_bytes = sizeof(unsigned int) * levels;

  constexpr static size_t content_memory_size_in_bytes =
      enumeration_x_size_in_bytes + coefficient_size_in_bytes + center_partsum_size_in_bytes +
      partdist_size_in_bytes + parent_indices_size_in_bytes + open_node_count_size_in_bytes;

public:
  __device__ __host__ inline SubtreeEnumerationBuffer(unsigned char *memory)
      : center_partsum(reinterpret_cast<enumf *>(memory)),
        partdist(reinterpret_cast<enumf *>(memory + center_partsum_size_in_bytes)),
        enumeration_x(reinterpret_cast<enumi *>(memory + center_partsum_size_in_bytes +
                                                partdist_size_in_bytes)),
        coefficients(reinterpret_cast<enumi *>(memory + center_partsum_size_in_bytes +
                                               partdist_size_in_bytes +
                                               enumeration_x_size_in_bytes)),
        parent_indices(reinterpret_cast<unsigned int *>(
            memory + center_partsum_size_in_bytes + partdist_size_in_bytes +
            enumeration_x_size_in_bytes + coefficient_size_in_bytes)),
        open_node_count(reinterpret_cast<unsigned int *>(
            memory + center_partsum_size_in_bytes + partdist_size_in_bytes +
            enumeration_x_size_in_bytes + coefficient_size_in_bytes + parent_indices_size_in_bytes))
  {
    assert(((intptr_t)memory) % sizeof(enumf) == 0);
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

  // ensure alignment
  constexpr static size_t memory_size_in_bytes =
      ((content_memory_size_in_bytes - 1) / sizeof(enumf) + 1) * sizeof(enumf);

  static_assert(memory_size_in_bytes >= content_memory_size_in_bytes,
                "Bug in memory_size_in_bytes calculation");

  static_assert(memory_size_in_bytes < std::numeric_limits<unsigned int>::max(),
                "Requires more memory than indexable with unsigned int");

  __device__ __host__ inline CudaEnumeration<dimensions_per_level>
  get_enumeration(unsigned int tree_level, unsigned int index, Matrix mu_block, const enumf *rdiag,
                  const uint32_t *radius_squared_location)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    CudaEnumeration<dimensions_per_level> result;

    const unsigned int offset_kk = (levels - tree_level - 1) * dimensions_per_level;

    result.mu                      = mu_block;
    result.rdiag                   = &rdiag[offset_kk];
    result.radius_squared_location = radius_squared_location;

    for (unsigned int i = 0; i < dimensions_per_level; ++i)
    {
      result.x[i] = enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                                  i * max_nodes_per_level + index];

      const enumf center_partsum_i = center_partsum[tree_level * dimensions * max_nodes_per_level +
                                                    (offset_kk + i) * max_nodes_per_level + index];
      assert(!isnan(center_partsum_i));
      result.center_partsums[i][dimensions_per_level - 1] = center_partsum_i;
    }

    result.center[dimensions_per_level - 1] =
        center_partsum[tree_level * dimensions * max_nodes_per_level +
                       (offset_kk + dimensions_per_level - 1) * max_nodes_per_level + index];
    result.partdist[dimensions_per_level - 1] = partdist[tree_level * max_nodes_per_level + index];
    return result;
  }

  __device__ __host__ inline void
  set_enumeration(unsigned int tree_level, unsigned int index,
                  const CudaEnumeration<dimensions_per_level> &value)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
  #ifndef NDEBUG
    const unsigned int offset_kk = (levels - tree_level - 1) * dimensions_per_level;
  #endif

    for (unsigned int i = 0; i < dimensions_per_level; ++i)
    {
      enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                    i * max_nodes_per_level + index] = value.x[i];

      assert(center_partsum[tree_level * dimensions * max_nodes_per_level + 
                            (offset_kk + i) * max_nodes_per_level + index] == 
             value.center_partsums[i][dimensions_per_level - 1]);
    }
    assert(center_partsum[tree_level * dimensions * max_nodes_per_level +
                          (offset_kk + dimensions_per_level - 1) * max_nodes_per_level + index] ==
           value.center[dimensions_per_level - 1]);
    assert(partdist[tree_level * max_nodes_per_level + index] ==
           value.partdist[dimensions_per_level - 1]);
  }

  __device__ __host__ inline void init_subtree(unsigned int tree_level, unsigned int index,
                                               enumf parent_partdist, enumf center)
  {
    for (unsigned int i = 0; i < dimensions_per_level; ++i)
    {
      enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                    i * max_nodes_per_level + index] = NAN;
    }
    partdist[tree_level * max_nodes_per_level + index] = parent_partdist;
    enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                  (dimensions_per_level - 1) * max_nodes_per_level + index] =
        static_cast<enumi>(round(center));
    assert(!isnan(static_cast<enumi>(round(center))));
  }

  __device__ __host__ inline void set_center_partsum(unsigned int tree_level, unsigned int index,
                                                     unsigned int orth_basis_index, enumf value)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    assert(orth_basis_index < dimensions);
    center_partsum[tree_level * dimensions * max_nodes_per_level +
                   orth_basis_index * max_nodes_per_level + index] = value;
  }

  __device__ __host__ inline enumf get_center_partsum(unsigned int tree_level, unsigned int index,
                                                      unsigned int orth_basis_index)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    assert(orth_basis_index < dimensions);
    return center_partsum[tree_level * dimensions * max_nodes_per_level +
                          orth_basis_index * max_nodes_per_level + index];
  }

  __device__ __host__ inline enumi get_coefficient(unsigned int tree_level, unsigned int index,
                                                   unsigned int coordinate)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    assert(coordinate < dimensions_per_level);
    return coefficients[tree_level * dimensions_per_level * max_nodes_per_level +
                        coordinate * max_nodes_per_level + index];
  }

  __device__ __host__ inline void set_coefficient(unsigned int tree_level, unsigned int index,
                                                  unsigned int coordinate, enumi value)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    assert(coordinate < dimensions_per_level);
    coefficients[tree_level * dimensions_per_level * max_nodes_per_level +
                 coordinate * max_nodes_per_level + index] = value;
  }

  __device__ __host__ inline unsigned int get_parent_index(unsigned int tree_level,
                                                           unsigned int index)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    return parent_indices[tree_level * max_nodes_per_level + index];
  }

  __device__ __host__ inline enumf get_partdist(unsigned int tree_level, unsigned int index)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    return partdist[tree_level * max_nodes_per_level + index];
  }

  __device__ __host__ inline void set_partdist(unsigned int tree_level, unsigned int index,
                                               enumf value)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    partdist[tree_level * max_nodes_per_level + index] = value;
  }

  __device__ __host__ inline unsigned int get_node_count(unsigned int tree_level)
  {
    assert(tree_level < levels);
    return open_node_count[tree_level];
  }

  __device__ __host__ inline unsigned int add_node(unsigned int tree_level,
                                                   unsigned int parent_node_index)
  {
    assert(tree_level < levels);
    const unsigned int new_task_index = aggregated_atomic_inc(&open_node_count[tree_level]);
    assert(new_task_index < max_nodes_per_level);
    // in this case, we want an error also in Release builds
    if (new_task_index >= max_nodes_per_level)
    {
      runtime_error();
    }
    parent_indices[tree_level * max_nodes_per_level + new_task_index] = parent_node_index;
    return new_task_index;
  }

  /**
   * Removes all nodes this functions was called for with keep_this_thread_task=false from the tree.
   *
   * To allow an efficient implementation, requires that old_index == node_count -
   * active_thread_count + cooperative_group.thread_rank(), i.e. all threads in the group have to
   * call this function for the last active_thread_count nodes in the tree. Calls with
   * cooperative_group.thread_rank() >= active_thread_count will be ignored.
   */
  template <typename CG, unsigned int block_size>
  __device__ __host__ inline void
  filter_nodes(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter,
               unsigned int tree_level, unsigned int old_index, bool keep_this_thread_task,
               unsigned int active_thread_count)
  {
    assert(tree_level < levels);
    assert(active_thread_count <= open_node_count[tree_level]);
    assert(old_index ==
           open_node_count[tree_level] - active_thread_count + cooperative_group.thread_rank());
    assert(tree_level + 1 == levels || open_node_count[tree_level + 1] == 0);

    unsigned int kept_tasks = 0;
    const bool is_active =
        keep_this_thread_task && cooperative_group.thread_rank() < active_thread_count;
    const unsigned int new_offset =
        prefix_counter.prefix_count(cooperative_group, is_active, kept_tasks);
    const unsigned int new_index = new_offset + open_node_count[tree_level] - active_thread_count;

    enumi coefficients_tmp[dimensions_per_level];
    enumf center_partsum_tmp[dimensions];
    enumf partdist_tmp;
    unsigned int parent_index_tmp;
    if (is_active)
    {
      partdist_tmp     = partdist[tree_level * max_nodes_per_level + old_index];
      parent_index_tmp = parent_indices[tree_level * max_nodes_per_level + old_index];
      for (unsigned int i = 0; i < dimensions_per_level; ++i)
      {
        coefficients_tmp[i] =
            enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                          i * max_nodes_per_level + old_index];
      }
      for (unsigned int i = 0; i < dimensions; ++i)
      {
        center_partsum_tmp[i] = center_partsum[tree_level * dimensions * max_nodes_per_level +
                                               i * max_nodes_per_level + old_index];
      }
    }

    cooperative_group.sync();

    if (is_active)
    {
      partdist[tree_level * max_nodes_per_level + new_index]       = partdist_tmp;
      parent_indices[tree_level * max_nodes_per_level + new_index] = parent_index_tmp;
      for (unsigned int i = 0; i < dimensions_per_level; ++i)
      {
        enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                      i * max_nodes_per_level + new_index] = coefficients_tmp[i];
      }
      for (unsigned int i = 0; i < dimensions; ++i)
      {
        center_partsum[tree_level * dimensions * max_nodes_per_level + i * max_nodes_per_level +
                       new_index] = center_partsum_tmp[i];
      }
    }

    if (cooperative_group.thread_rank() == 0)
    {
      open_node_count[tree_level] -= active_thread_count - kept_tasks;
    }
  }
};

template <typename eval_sol_fn, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
struct ProcessLeafCallback
{
  unsigned int level;
  unsigned int parent_index;
  unsigned int start_point_dim;
  eval_sol_fn &process_sol;
  Matrix mu;
  const enumi *start_points;
  uint32_t *radius_squared_location;
  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer;

  __device__ __host__ void operator()(const enumi *x, enumf squared_norm);
};

template <typename eval_sol_fn, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ inline void
ProcessLeafCallback<eval_sol_fn, levels, dimensions_per_level, max_nodes_per_level>::operator()(
    const enumi *x, enumf squared_norm)
{
  if (squared_norm == 0)
  {
    return;
  }
  for (unsigned int i = 0; i < dimensions_per_level; ++i)
  {
    process_sol(x[i], i, false, squared_norm, radius_squared_location);
  }
  unsigned int index = parent_index;
  for (unsigned int j = levels - 1; j > 0; --j)
  {
    for (unsigned int i = 0; i < dimensions_per_level; ++i)
    {
      process_sol(buffer.get_coefficient(j, index, i), i + (levels - j) * dimensions_per_level,
                  false, squared_norm, radius_squared_location);
    }
    index = buffer.get_parent_index(j, index);
  }
  index = buffer.get_parent_index(0, index);
  for (unsigned int i = 0; i < start_point_dim; ++i)
  {
    process_sol(start_points[index * start_point_dim + i], i + levels * dimensions_per_level,
                i + 1 == start_point_dim, squared_norm, radius_squared_location);
  }
}

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct AddToTreeCallback
{
  unsigned int level;
  unsigned int parent_index;
  Matrix mu;
  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer;
  PerfCounter &counter;

  __device__ __host__ void operator()(const enumi *x, enumf squared_norm);
};

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void
AddToTreeCallback<levels, dimensions_per_level, max_nodes_per_level>::operator()(const enumi *x,
                                                                                 enumf squared_norm)
{
  assert(level > 0);

  const unsigned int new_index = buffer.add_node(level, parent_index);
  for (unsigned int j = 0; j < dimensions_per_level; ++j)
  {
    buffer.set_coefficient(level, new_index, j, x[j]);
  }
  buffer.set_partdist(level, new_index, squared_norm);
  // subtree initialization will be done later in a synchronized way
}

/**
 * Calculates the difference of this center partsum to the center partsum of the parent point
 */
template <unsigned int levels, unsigned int dimensions_per_level>
__device__ __host__ inline enumf calc_center_partsum_delta(unsigned int level, unsigned int index,
                                                           unsigned int center_partsum_index,
                                                           enumi x[dimensions_per_level], Matrix mu)
{
  unsigned int kk_offset = (levels - level - 1) * dimensions_per_level;
  enumf center_partsum = 0;
  for (unsigned int j = 0; j < dimensions_per_level; ++j)
  {
    center_partsum -= x[j] * mu.at(center_partsum_index, j + dimensions_per_level + kk_offset);
  }
  assert(!isnan(center_partsum));
  return center_partsum;
}

/**
 * Initializes newly generated nodes with all information necessary to perform subtree enumeration,
 * namely the center_partsums and the partdist. Requires the coefficients of these nodes to be
 * already stored in the tree buffer.
 *
 * Newly generated nodes are all nodes on the given level, except the already_calculated_node_count
 * first nodes.
 */
template <typename CG, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ inline void
init_new_nodes(CG &group, unsigned int level, unsigned int already_calculated_node_count,
               SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
               Matrix mu, PerfCounter &counter)
{
  for (unsigned int new_index = already_calculated_node_count + group.thread_rank();
       new_index < buffer.get_node_count(level); new_index += group.size())
  {
    unsigned int kk_offset = (levels - level - 1) * dimensions_per_level;
    unsigned int center_i  = kk_offset + dimensions_per_level - 1;

    const unsigned int parent_index = buffer.get_parent_index(level, new_index);
    enumi x[dimensions_per_level];
    for (unsigned int j = 0; j < dimensions_per_level; ++j)
    {
      x[j] = buffer.get_coefficient(level, new_index, j);
    }

    // sets center_partsum[i] = parent_center_partsum[i] + calc_center_partsum_delta(..., i)
    // to reduce latency, the loop is transformed as to load data now that is needed after some loop
    // cycles
    constexpr unsigned int loop_preload_count  = 3;
    constexpr unsigned int loop_preload_offset = loop_preload_count - 1;
    unsigned int i                             = 0;
    enumf center_partsum;
    enumf preloaded_parent_center_partsums[loop_preload_count];

#pragma unroll
    for (unsigned int j = 0; j < loop_preload_offset; ++j)
    {
      preloaded_parent_center_partsums[j] = buffer.get_center_partsum(level - 1, parent_index, j);
    }

    for (; i + 2 * loop_preload_offset < kk_offset + dimensions_per_level; i += loop_preload_count)
    {
#pragma unroll
      for (unsigned int j = 0; j < loop_preload_count; ++j)
      {
        preloaded_parent_center_partsums[(j + loop_preload_offset) % loop_preload_count] =
            buffer.get_center_partsum(level - 1, parent_index, i + j + loop_preload_offset);

        assert(preloaded_parent_center_partsums[j] ==
               buffer.get_center_partsum(level - 1, parent_index, i + j));
        center_partsum =
            preloaded_parent_center_partsums[j] +
            calc_center_partsum_delta<levels, dimensions_per_level>(level, new_index, i + j, x, mu);
        buffer.set_center_partsum(level, new_index, i + j, center_partsum);
      }
    }

    assert(i + 2 * loop_preload_offset - loop_preload_count + 1 <=
           kk_offset + dimensions_per_level);
    assert(i + 2 * loop_preload_offset >= kk_offset + dimensions_per_level);

#pragma unroll
    for (unsigned int ii = 2 * loop_preload_offset - loop_preload_count + 1;
         ii <= 2 * loop_preload_offset; ++ii)
    {
      if (i + ii == kk_offset + dimensions_per_level)
      {
#pragma unroll
        for (unsigned int j = 0; j < ii; ++j)
        {
          if (j + loop_preload_offset < ii)
          {
            preloaded_parent_center_partsums[(j + loop_preload_offset) % loop_preload_count] =
                buffer.get_center_partsum(level - 1, parent_index, i + j + loop_preload_offset);
          }
          assert(preloaded_parent_center_partsums[j % loop_preload_count] ==
                 buffer.get_center_partsum(level - 1, parent_index, i + j));
          center_partsum = preloaded_parent_center_partsums[j % loop_preload_count] +
                           calc_center_partsum_delta<levels, dimensions_per_level>(level, new_index,
                                                                                   i + j, x, mu);
          buffer.set_center_partsum(level, new_index, i + j, center_partsum);
        }
      }
    }

    enumf center = buffer.get_center_partsum(level, new_index, center_i);
    assert(!isnan(center));
    buffer.set_center_partsum(level, new_index, center_i, center);

    enumf partdist = buffer.get_partdist(level, new_index);
    buffer.init_subtree(level, new_index, partdist, center);
  }
}

/**
 * Generates more children for the last group.size() nodes on the given level and adds them to the
 * buffer. The subtree enumerations of the processed nodes are accordingly updated so that they will
 * only yield new children.
 */
template <typename CG, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ void generate_nodes_children(
    CG &group, SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
    int level, Matrix mu, const enumf *rdiag, uint32_t *radius_squared_location,
    unsigned int max_subtree_paths, PerfCounter &counter)
{
  assert(level < levels - 1);
  assert(level >= 0);

  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  unsigned int max_paths       = max_subtree_paths;
  const unsigned int offset_kk = (levels - level - 1) * dimensions_per_level;
  unsigned int existing_nodes  = buffer.get_node_count(level + 1);

  group.sync();

  if (active)
  {
    CudaEnumeration<dimensions_per_level> enumeration = buffer.get_enumeration(
        level, index, mu.block(offset_kk, offset_kk), rdiag, radius_squared_location);

    bool is_done = enumeration.template is_enumeration_done<dimensions_per_level - 1>();

    if (!is_done)
    {
      typedef AddToTreeCallback<levels, dimensions_per_level, max_nodes_per_level> CallbackType;
      CallbackType callback = {static_cast<unsigned int>(level + 1), index, mu, buffer, counter};
      enumeration.template enumerate_recursive(
          callback, max_paths, counter, kk_marker<dimensions_per_level - 1>());

      buffer.set_enumeration(level, index, enumeration);
    }

  }

  group.sync();

  init_new_nodes(group, level + 1, existing_nodes, buffer, mu, counter);
}

/**
 * Searches the subtrees of the last group.size() nodes on the last tree level, possibly finding a
 * new nonzero shortest vector. The subtree enumerations of the processed nodes are accordingly
 * updated so that they will only yield new vectors.
 */
template <typename CG, typename eval_sol_fn, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ void inline process_leaf_nodes(
    CG &group, SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
    Matrix mu, const enumf *rdiag, uint32_t *radius_squared_location, unsigned int max_paths,
    eval_sol_fn &process_sol, const enumi *start_points, unsigned int start_point_dim,
    PerfCounter &node_counter)
{
  const unsigned int level               = levels - 1;
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  if (active)
  {
    CudaEnumeration<dimensions_per_level> enumeration =
        buffer.get_enumeration(level, index, mu, rdiag, radius_squared_location);

    typedef ProcessLeafCallback<eval_sol_fn, levels, dimensions_per_level, max_nodes_per_level>
        CallbackT;
    CallbackT callback = {level + 1, index,        start_point_dim,         process_sol,
                          mu,        start_points, radius_squared_location, buffer};
    enumeration.template enumerate_recursive(
        callback, max_paths, node_counter, kk_marker<dimensions_per_level - 1>());

    buffer.set_enumeration(level, index, enumeration);
  }
}

/**
 * Calculates the count of nodes among the last roup.size() nodes on the given level whose subtrees
 * have nodes not exceeding the radius limit.
 */
template <typename CG, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ inline unsigned int get_done_node_count(
    CG &group, unsigned int *shared_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer, int level,
    Matrix mu, const enumf *rdiag, const uint32_t *radius_square_location)
{
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  const unsigned int offset_kk = (levels - level - 1) * dimensions_per_level;

  *shared_counter = 0;

  group.sync();

  if (active)
  {
    bool is_done = buffer
                       .get_enumeration(level, index, mu.block(offset_kk, offset_kk), rdiag,
                                        radius_square_location)
                       .template is_enumeration_done<dimensions_per_level - 1>();
    if (is_done)
    {
      aggregated_atomic_inc(shared_counter);
    }
  }

  group.sync();

  return *shared_counter;
}

/**
 * Removes all nodes among the last group.size() nodes on the given level whose subtrees have nodes
 * with partdist exceeding the radius limit. Be careful as the buffer can still have nodes
 * referencing such a done node as a parent node, since the enumeration data is updated when
 * children are generated, not when children are fully processed.
 */
template <typename CG, unsigned int block_size, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void remove_done_nodes(
    CG &group, PrefixCounter<CG, block_size> &prefix_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer, int level,
    Matrix mu, const enumf *rdiag, const uint32_t *radius_square_location)
{
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  const unsigned int offset_kk = (levels - level - 1) * dimensions_per_level;

  bool is_done = buffer
                     .get_enumeration(level, index, mu.block(offset_kk, offset_kk), rdiag,
                                      radius_square_location)
                     .template is_enumeration_done<dimensions_per_level - 1>();

  group.sync();

  buffer.filter_nodes(group, prefix_counter, level, index, !is_done, active_thread_count);
}

struct StrategyOpts
{
  unsigned int max_subtree_paths;
  // stop children generation when the percentage of parent points that have still unprocessed
  // children drops beneath this percentage
  float min_active_parents_percentage;
  // stop children generation when the count of children exceeds this limit
  unsigned int max_children_buffer_size;
};

/**
 * Generates children from the last group.size() nodes on level and adds them to the buffer, until
 * either the children buffer is full or most of these nodes are done.
 */
template <typename CG, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ inline void generate_level_children(
    CG &group, unsigned int *shared_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer, int level,
    Matrix mu, const enumf *rdiag, uint32_t *radius_square_location, PerfCounter counter,
    StrategyOpts opts)
{
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());

  while (true)
  {
    generate_nodes_children(group, buffer, level, mu, rdiag, radius_square_location,
                            opts.max_subtree_paths, counter);

    group.sync();

    const unsigned int done_node_count = get_done_node_count(group, shared_counter, buffer, level,
                                                             mu, rdiag, radius_square_location);

    group.sync();

    if (CUENUM_TRACE && thread_id() == 0)
    {
      printf("Thread 0: Worked on level %d, next level points are %d, %d nodes of current working pool (%d) are done\n",
             level, buffer.get_node_count(level + 1), *shared_counter, active_thread_count);
    }

    if (buffer.get_node_count(level + 1) >= opts.max_children_buffer_size)
    {
      break;
    }
    else if (done_node_count >= active_thread_count * (1 - opts.min_active_parents_percentage))
    {
      break;
    }
    group.sync();
  }
}

template <typename CG, unsigned int block_size, typename eval_sol_fn, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void process_leaf_level(
    CG &group, PrefixCounter<CG, block_size> &prefix_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer, Matrix mu,
    const enumf *rdiag, uint32_t *radius_square_location, PerfCounter node_counter,
    eval_sol_fn &process_sol, const enumi *start_points, unsigned int start_point_dim,
    StrategyOpts opts)
{
  const unsigned int level = levels - 1;
  while (buffer.get_node_count(level) > 0)
  {
    process_leaf_nodes(group, buffer, mu, rdiag, radius_square_location, opts.max_subtree_paths,
                       process_sol, start_points, start_point_dim, node_counter);
    remove_done_nodes(group, prefix_counter, buffer, level, mu, rdiag, radius_square_location);
    group.sync();
  }
}

/**
 * Removes finished nodes from the parent level of the given level. Does nothing when called on the
 * root.
 */
template <typename CG, unsigned int block_size, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void cleanup_parent_level(
    CG &group, PrefixCounter<CG, block_size> &prefix_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer, int level,
    Matrix mu, const enumf *rdiag, uint32_t *radius_square_location)
{
  if (level > 0)
  {
    group.sync();

    remove_done_nodes(group, prefix_counter, buffer, level - 1, mu, rdiag, radius_square_location);

    group.sync();

    if (CUENUM_TRACE && thread_id() == 0)
    {
      printf("Thread 0: Cleaned up level %d, has now %d nodes\n", level - 1, buffer.get_node_count(level - 1));
    }

    group.sync();
  }
}

template <typename CG, unsigned int block_size, typename eval_sol_fn, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void
clear_level(CG &group, PrefixCounter<CG, block_size> &prefix_counter, unsigned int *shared_counter,
            SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
            int level, Matrix mu, const enumf *rdiag, uint32_t *radius_square_location,
            eval_sol_fn &process_sol, const enumi *start_points, unsigned int start_point_dim,
            StrategyOpts opts, PerfCounter node_counter)
{
  while (level >= 0)
  {
    if (level + 1 < levels)
    {
      if (buffer.get_node_count(level) > 0)
      {
        generate_level_children(group, shared_counter, buffer, level, mu, rdiag,
                                radius_square_location, node_counter, opts);
        ++level;
      }
      else
      {
        cleanup_parent_level(group, prefix_counter, buffer, level, mu, rdiag,
                             radius_square_location);
        --level;
      }
    }
    else
    {
      process_leaf_level(group, prefix_counter, buffer, mu, rdiag, radius_square_location,
                         node_counter, process_sol, start_points, start_point_dim, opts);
      cleanup_parent_level(group, prefix_counter, buffer, level, mu, rdiag, radius_square_location);
      --level;
    }
    group.sync();
  }
}

constexpr unsigned int enumerate_block_size               = 128;
constexpr unsigned int enumerate_cooperative_group_size   = 32;
constexpr unsigned int enumerate_point_stream_buffer_size = 100;

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct Opts
{
  StrategyOpts tree_clear_opts;
  unsigned int initial_nodes_per_group;
  unsigned int thread_count;
};

template <unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__global__ void __launch_bounds__(enumerate_block_size, 2)
    enumerate_kernel(unsigned char *buffer_memory, const enumi *start_points,
                     unsigned int *processed_start_point_counter, unsigned int start_point_count,
                     unsigned int start_point_dim, const enumf *mu_ptr, const enumf *rdiag,
                     uint32_t *radius_squared_location, unsigned long long *perf_counter_memory,
                     unsigned char* point_stream_memory,
                     Opts<levels, dimensions_per_level, max_nodes_per_level> opts)
{
  typedef cooperative_groups::thread_block_tile<32> CG;
  typedef SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> SubtreeBuffer;
  typedef PointStreamEvaluator<enumerate_point_stream_buffer_size> Evaluator;

  constexpr unsigned int block_size            = enumerate_block_size;
  constexpr unsigned int dimensions            = dimensions_per_level * levels;
  constexpr unsigned int group_count_per_block = enumerate_block_size / enumerate_cooperative_group_size;

  constexpr unsigned int mu_shared_memory_size    = dimensions * dimensions * sizeof(enumf);
  constexpr unsigned int rdiag_shared_memory_size = dimensions * sizeof(enumf);

  constexpr unsigned int shared_mem_size = group_count_per_block * sizeof(unsigned int) + sizeof(unsigned int) + 
                                           mu_shared_memory_size + rdiag_shared_memory_size;

  __shared__ unsigned char shared_mem[shared_mem_size];

  CG group = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());
  const unsigned int group_id          = thread_id() / enumerate_cooperative_group_size;
  const unsigned int group_id_in_block = thread_id_in_block() / enumerate_cooperative_group_size;
  assert(group.size() == enumerate_cooperative_group_size);

  PrefixCounter<CG, block_size> prefix_counter;

  enumf *mu_shared    = reinterpret_cast<enumf *>(shared_mem);
  enumf *rdiag_shared = reinterpret_cast<enumf *>(shared_mem + mu_shared_memory_size);

  unsigned int *group_shared_counter =
      reinterpret_cast<unsigned int *>(shared_mem + group_id_in_block * sizeof(unsigned int) +
                                       mu_shared_memory_size + rdiag_shared_memory_size);
  unsigned int *point_stream_counter = 
      reinterpret_cast<unsigned int *>(shared_mem + group_count_per_block * sizeof(unsigned int) +
                                       mu_shared_memory_size + rdiag_shared_memory_size);

  const unsigned int ldmu = dimensions + start_point_dim;
  for (unsigned int i = threadIdx.x; i < dimensions * dimensions; i += blockDim.x)
  {
    mu_shared[i] = mu_ptr[i / dimensions * ldmu + i % dimensions];
  }
  for (unsigned int i = threadIdx.x; i < dimensions; i += blockDim.x)
  {
    rdiag_shared[i] = rdiag[i];
  }
  __syncthreads();
  Matrix mu(mu_shared, dimensions);

  PerfCounter node_counter(perf_counter_memory);

  Evaluator process_sol(point_stream_memory + blockIdx.x * Evaluator::memory_size_in_bytes(dimensions + start_point_dim), point_stream_counter);

  assert(opts.initial_nodes_per_group <= group.size());

  SubtreeBuffer buffer(buffer_memory + group_id * SubtreeBuffer::memory_size_in_bytes);

  while (true)
  {
    group.sync();
    if (group.thread_rank() == 0)
    {
      *group_shared_counter =
          atomic_add(processed_start_point_counter, opts.initial_nodes_per_group);
    }
    buffer.init(group);
    group.sync();

    if (*group_shared_counter >= start_point_count)
    {
      break;
    }
    const unsigned int start_point_index = *group_shared_counter + group.thread_rank();
    const bool active =
        group.thread_rank() < opts.initial_nodes_per_group && start_point_index < start_point_count;
    if (active)
    {
      const enumi *start_point = &start_points[start_point_index * start_point_dim];
      const unsigned int index = buffer.add_node(0, start_point_index);
      for (unsigned int i = 0; i < dimensions; ++i)
      {
        enumf center_partsum = 0;
        for (unsigned int j = 0; j < start_point_dim; ++j)
        {
          center_partsum -= start_point[j] * mu_ptr[i * ldmu + dimensions + j];
        }
        assert(!isnan(center_partsum));
        buffer.set_center_partsum(0, index, i, center_partsum);
      }
      enumf partdist = 0;
      for (int j = 0; j < start_point_dim; ++j)
      {
        enumf alpha = start_point[j];
        for (unsigned int i = j + 1; i < start_point_dim; ++i)
        {
          alpha += start_point[i] * mu_ptr[(j + dimensions) * ldmu + i + dimensions];
        }
        assert(rdiag[dimensions + j] >= 0);
        partdist += alpha * alpha * rdiag[dimensions + j];
        assert(partdist >= 0);
      }
      buffer.init_subtree(0, index, partdist, buffer.get_center_partsum(0, index, dimensions - 1));
    }
    if (CUENUM_TRACE && thread_id() == 0)
    {
      printf("Thread 0: Get %d new nodes\n", opts.initial_nodes_per_group);
    }

    group.sync();

    clear_level(group, prefix_counter, group_shared_counter, buffer, 0, mu, rdiag_shared,
                radius_squared_location, process_sol, start_points, start_point_dim,
                opts.tree_clear_opts, node_counter);
  }
}

constexpr unsigned int get_grid_size(unsigned int thread_count)
{
  return (thread_count - 1) / enumerate_block_size + 1;
}

constexpr unsigned int get_started_thread_count(unsigned int thread_count) {
  return get_grid_size(thread_count) * enumerate_block_size;
}
 
typedef std::function<float(double, float*)> process_sol_fn;

/**
 * Enumerates all points within the enumeration bound (initialized to parameter initial_radius) and calls
 * the given function on each coordinate of each of these points.
 * 
 * For the description of the parameters, have the total lattice dimension n = levels * dimensions_per_level + start_point_dim.
 * 
 * @param mu - pointer to memory containing the normalized gram-schmidt coefficients, in row-major format.
 * In other words, the memory must contain n consecutive batches of memory, each consisting of n entries storing
 * the values of the corresponding row of the matrix. 
 * @param rdiag - n entries containing the squared norms of the gram-schmidt vectors, in one contigous segment of memory
 * @param start_points - Function yielding a pointer to memory containing the start_point_dim coefficients of the i-th start point. This
 * pointer must stay valid until the next call of the function. 
 * @param process_sol - callback function called on solution points that were found. This is a host function.
 */
template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level, bool print_status = true>
void enumerate(const enumf *mu, const enumf *rdiag, const float *start_points,
               unsigned int start_point_dim, unsigned int start_point_count, enumf initial_radius,
               process_sol_fn process_sol,
               Opts<levels, dimensions_per_level, max_nodes_per_level> opts)
{
  typedef SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> SubtreeBuffer;
  typedef PointStreamEvaluator<enumerate_point_stream_buffer_size> Evaluator;
  typedef PointStreamEndpoint<enumerate_point_stream_buffer_size> PointStream;

  constexpr unsigned int tree_dimensions = levels * dimensions_per_level;
  const unsigned int mu_n                = tree_dimensions + start_point_dim;
  const unsigned int grid_size           = get_grid_size(opts.thread_count);
  const unsigned int group_count         = grid_size * enumerate_block_size / enumerate_cooperative_group_size;

  CudaPtr<unsigned char> buffer_mem =
      cuda_alloc(unsigned char, SubtreeBuffer::memory_size_in_bytes * group_count);
  CudaPtr<unsigned char> point_stream_memory = 
      cuda_alloc(unsigned char, Evaluator::memory_size_in_bytes(mu_n) * grid_size);

  CudaPtr<uint32_t> radius_mem             = cuda_alloc(uint32_t, 1);
  CudaPtr<enumf> device_mu                 = cuda_alloc(enumf, mu_n * mu_n);
  CudaPtr<enumf> device_rdiag              = cuda_alloc(enumf, mu_n);
  CudaPtr<unsigned long long> node_counter = cuda_alloc(unsigned long long, 1);
  CudaPtr<enumi> device_start_points       = cuda_alloc(enumi, start_point_count * start_point_dim);
  CudaPtr<unsigned int> processed_start_point_count = cuda_alloc(unsigned int, 1);

  const uint32_t repr_initial_radius_squared =
      float_to_int_order_preserving_bijection(initial_radius * initial_radius);
  check(cudaMemcpy(device_mu.get(), mu, mu_n * mu_n * sizeof(enumf), cudaMemcpyHostToDevice));
  check(cudaMemcpy(device_rdiag.get(), rdiag, mu_n * sizeof(enumf), cudaMemcpyHostToDevice));
  check(cudaMemcpy(radius_mem.get(), &repr_initial_radius_squared, sizeof(uint32_t),
                   cudaMemcpyHostToDevice));
  check(cudaMemcpy(device_start_points.get(), start_points,
                   start_point_dim * start_point_count * sizeof(enumi), cudaMemcpyHostToDevice));

  PointStream stream(point_stream_memory.get(), radius_mem.get(), grid_size, mu_n);
  stream.init();

  cudaEvent_t raw_event;
  check(cudaEventCreateWithFlags(&raw_event, cudaEventDisableTiming));
  CudaEvent event(raw_event);

  cudaStream_t raw_exec_stream;
  check(cudaStreamCreate(&raw_exec_stream));
  CudaStream exec_stream(raw_exec_stream);

  if (print_status)
  {
    std::cout << "Enumerating " << (levels * dimensions_per_level)
              << " dimensional lattice using cuda, started " << grid_size << " block with "
              << enumerate_block_size << " threads each" << std::endl;
  }
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

  enumerate_kernel<<<dim3(grid_size), dim3(enumerate_block_size), 0, exec_stream.get()>>>(
      buffer_mem.get(), device_start_points.get(), processed_start_point_count.get(),
      start_point_count, start_point_dim, device_mu.get(), device_rdiag.get(), radius_mem.get(),
      node_counter.get(), point_stream_memory.get(), opts);
  
  check(cudaEventRecord(event.get(), exec_stream.get()));

  while(cudaEventQuery(event.get()) != cudaSuccess) 
  {
    stream.query_new_points<process_sol_fn, print_status>(process_sol);
  }
  stream.wait_for_event(event.get());
  stream.query_new_points<process_sol_fn, print_status>(process_sol);

  check(cudaDeviceSynchronize());
  check(cudaGetLastError());

  if (print_status)
  {
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Enumeration done in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;

    unsigned long long searched_nodes;
    uint32_t result_radius;
    check(cudaMemcpy(&searched_nodes, node_counter.get(), sizeof(unsigned long long),
                     cudaMemcpyDeviceToHost));

    check(cudaMemcpy(&result_radius, radius_mem.get(), sizeof(uint32_t), cudaMemcpyDeviceToHost));
    std::cout << "Searched " << searched_nodes
              << " tree nodes, and decreased enumeration bound down to "
              << sqrt(int_to_float_order_preserving_bijection(result_radius)) << std::endl;
  }
}

}

#endif