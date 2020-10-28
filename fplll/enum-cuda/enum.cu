#include "cooperative_groups.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <array>
#include <assert.h>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <stdint.h>
#include <vector>

#include "atomic.cuh"
#include "prefix.cuh"
#include "recenum.cuh"
#include "util.h"

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct SubtreeEnumerationBuffer
{
private:
  // shape [levels, dimensions_per_level, max_nodes_per_level]
  enumi *enumeration_x;
  // shape [levels, dimensions_per_level, max_nodes_per_level]
  enumi *coefficients;
  // shape [levels, dimensions, max_nodes_per_level], of subtree root
  enumf *center_partsum;
  // shape [levels, max_nodes_per_level], of subtree root
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

  constexpr static size_t memory_size_in_bytes =
      enumeration_x_size_in_bytes + coefficient_size_in_bytes + center_partsum_size_in_bytes +
      partdist_size_in_bytes + parent_indices_size_in_bytes + open_node_count_size_in_bytes;

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
    const unsigned int offset_kk = (levels - tree_level - 1) * dimensions_per_level;

    for (unsigned int i = 0; i < dimensions_per_level; ++i)
    {
      enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                    i * max_nodes_per_level + index] = value.x[i];

      const enumf old_parent_center_partsum =
          center_partsum[tree_level * dimensions * max_nodes_per_level +
                         (offset_kk + i) * max_nodes_per_level + index];
      const enumf new_parent_center_partsum = value.center_partsums[i][dimensions_per_level - 1];
      assert(!isnan(old_parent_center_partsum));
      assert(!isnan(new_parent_center_partsum));
      assert(old_parent_center_partsum == new_parent_center_partsum);
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
    partdist[tree_level * max_nodes_per_level + index]                      = parent_partdist;
    enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                  (dimensions_per_level - 1) * max_nodes_per_level + index] = round(center);
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
#ifdef __CUDA_ARCH__
      asm("trap;");
#else
      throw;
#endif
    }
    parent_indices[tree_level * max_nodes_per_level + new_task_index] = parent_node_index;
    return new_task_index;
  }

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

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct ProcessLeafCallback
{
  unsigned int level;
  unsigned int parent_index;
  Matrix mu;
  uint32_t *radius_squared_location;
  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer;

  __device__ __host__ void operator()(const enumi *x, enumf squared_norm);
};

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void
ProcessLeafCallback<levels, dimensions_per_level, max_nodes_per_level>::operator()(
    const enumi *x, enumf squared_norm)
{
  if (squared_norm == 0)
  {
    return;
  }

  uint32_t squared_norm_repr = float_to_int_order_preserving_bijection(squared_norm);
  uint32_t old_repr          = atomic_min(radius_squared_location, squared_norm_repr);

  if (old_repr > squared_norm_repr)
  {
    // Here save the found result
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
  enumf center_partsum   = 0;
  for (unsigned int j = 0; j < dimensions_per_level; ++j)
  {
    center_partsum -= x[j] * mu.at(center_partsum_index, j + dimensions_per_level + kk_offset);
  }
  assert(!isnan(center_partsum));
  return center_partsum;
}

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
    // to reduce latency, the loop is transformed as to load data now that is needed after 4 loop cycles
    constexpr unsigned int loop_preload_count = 3;
    constexpr unsigned int loop_preload_offset = loop_preload_count - 1;
    unsigned int i = 0;
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

  unsigned long long end = time();
}

template <typename CG, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ void generate_nodes_children(
    CG &group, SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
    int level, Matrix mu, const enumf *rdiag, uint32_t *radius_squared_location,
    unsigned int max_subtree_paths, PerfCounter &counter)
{
  assert(level < levels - 1);

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

    typedef AddToTreeCallback<levels, dimensions_per_level, max_nodes_per_level> CallbackType;
    CallbackType callback = {level + 1, index, mu, buffer, counter};
    enumeration.template enumerate_recursive<dimensions_per_level - 1, CallbackType>(
        callback, max_paths, counter);

    buffer.set_enumeration(level, index, enumeration);
  }

  group.sync();

  init_new_nodes(group, level + 1, existing_nodes, buffer, mu, counter);
}

template <typename CG, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ void inline process_leaf_nodes(
    CG &group, SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
    Matrix mu, const enumf *rdiag, uint32_t *radius_squared_location, unsigned int max_paths,
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

    typedef ProcessLeafCallback<levels, dimensions_per_level, max_nodes_per_level> CallbackT;
    CallbackT callback = {level + 1, index, mu, radius_squared_location, buffer};
    enumeration.template enumerate_recursive<dimensions_per_level - 1, CallbackT>(
        callback, max_paths, node_counter);

    buffer.set_enumeration(level, index, enumeration);
  }
}

template <typename CG, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ inline unsigned int get_done_node_count(
    CG &group, unsigned int *shared_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
    int level, Matrix mu, const enumf *rdiag, const uint32_t *radius_square_location)
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

template <typename CG, unsigned int block_size, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void remove_done_nodes(
    CG &group, PrefixCounter<CG, block_size> &prefix_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
    int level, Matrix mu, const enumf *rdiag, const uint32_t *radius_square_location)
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
 * Generates children from the last group.size() nodes on level, until either
 * the children buffer is full or most of these nodes are done.
 */
template <typename CG, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
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

    debug_message_thread(
        "Worked on level %d, next level points are %d, %d nodes of current working pool "
        "(%d) are done\n",
        level, buffer.get_node_count(level + 1), *shared_counter, active_thread_count);

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

template <typename CG, unsigned int block_size, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void process_leaf_level(
    CG &group, PrefixCounter<CG, block_size> &prefix_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer, Matrix mu,
    const enumf *rdiag, uint32_t *radius_square_location, PerfCounter node_counter,
    StrategyOpts opts)
{
  const unsigned int level = levels - 1;
  while (buffer.get_node_count(level) > 0)
  {
    process_leaf_nodes(group, buffer, mu, rdiag, radius_square_location, opts.max_subtree_paths, node_counter);
    remove_done_nodes(group, prefix_counter, buffer, level, mu, rdiag, radius_square_location);
    group.sync();
  }
}

template <typename CG, unsigned int block_size, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void
cleanup_level(CG &group, PrefixCounter<CG, block_size> &prefix_counter,
              SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
              int level, Matrix mu, const enumf *rdiag, uint32_t *radius_square_location)
{
  if (level >= 0)
  {
    group.sync();

    remove_done_nodes(group, prefix_counter, buffer, level, mu, rdiag, radius_square_location);

    group.sync();

    debug_message_thread("Cleaned up level %d, has now %d nodes\n", level,
                         buffer.get_node_count(level));

    group.sync();
  }
}

template <typename CG, unsigned int block_size, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void
clear_level(CG &group, PrefixCounter<CG, block_size> &prefix_counter, unsigned int *shared_counter,
            SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
            int level, Matrix mu, const enumf *rdiag, uint32_t *radius_square_location,
            StrategyOpts opts, PerfCounter node_counter)
{
  while (level >= 0)
  {
    assert(all_threads_eq(group, level, shared_counter));
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
        --level;
        cleanup_level(group, prefix_counter, buffer, level, mu, rdiag, radius_square_location);
      }
    }
    else
    {
      process_leaf_level(group, prefix_counter, buffer, mu, rdiag, radius_square_location,
                         node_counter, opts);
      --level;
      cleanup_level(group, prefix_counter, buffer, level, mu, rdiag, radius_square_location);
    }
    group.sync();
  }
}

constexpr __device__ __host__ unsigned int constexpr_max(unsigned int a, unsigned int b)
{
  if (a > b)
  {
    return a;
  }
  else
  {
    return b;
  }
}

constexpr unsigned int search_block_size = 128;

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct Opts
{
  StrategyOpts tree_clear_opts;
  unsigned int initial_nodes_per_group;
  unsigned int start_point_dim;
  unsigned int thread_count;
};

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__global__ void __launch_bounds__(search_block_size, 1)
    search_kernel(unsigned char *buffer_memory, const enumi *start_points,
                  unsigned int *processed_start_point_counter, unsigned int start_point_count,
                  const enumf *mu_ptr, const enumf *rdiag, uint32_t *radius_squared_location,
                  unsigned long long *perf_counter_memory,
                  Opts<levels, dimensions_per_level, max_nodes_per_level> opts)
{
  typedef cooperative_groups::thread_block_tile<32> CG;
  typedef SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> SubtreeBuffer;

  constexpr unsigned int block_size            = search_block_size;
  constexpr unsigned int dimensions            = dimensions_per_level * levels;
  constexpr unsigned int group_count_per_block = block_size / 32;

  constexpr unsigned int mu_shared_memory_size = dimensions * dimensions * sizeof(enumf);

  constexpr unsigned int shared_mem_size =
      group_count_per_block * sizeof(unsigned int) + mu_shared_memory_size;

  extern __shared__ unsigned char shared_mem[shared_mem_size];

  CG group = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());
  const unsigned int group_id          = thread_id() / 32;
  const unsigned int group_id_in_block = thread_id_in_block() / 32;

  PrefixCounter<CG, block_size> prefix_counter;

  unsigned int *group_shared_counter =
      reinterpret_cast<unsigned int *>(shared_mem + group_id_in_block * sizeof(unsigned int));

  enumf *mu_shared =
      reinterpret_cast<enumf *>(shared_mem + group_count_per_block * sizeof(unsigned int));

  const unsigned int ldmu = dimensions + opts.start_point_dim;
  for (unsigned int i = threadIdx.x; i < dimensions * dimensions; i += blockDim.x)
  {
    mu_shared[i] = mu_ptr[i / dimensions * ldmu + i % dimensions];
  }
  __syncthreads();
  Matrix mu(mu_shared, dimensions);

  PerfCounter node_counter(perf_counter_memory);

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
      const enumi *start_point = &start_points[start_point_index * opts.start_point_dim];
      const unsigned int index = buffer.add_node(0, start_point_index);
      for (unsigned int i = 0; i < dimensions; ++i)
      {
        enumf center_partsum = 0;
        for (unsigned int j = 0; j < opts.start_point_dim; ++j)
        {
          center_partsum -= start_point[j] * mu_ptr[i * ldmu + dimensions + j];
        }
        buffer.set_center_partsum(0, index, i, center_partsum);
      }
      enumf partdist = 0;
      for (int j = 0; j < opts.start_point_dim; ++j)
      {
        enumf alpha = start_point[j];
        for (unsigned int i = j + 1; i < opts.start_point_dim; ++i)
        {
          alpha += start_point[i] * mu_ptr[(j + dimensions) * ldmu + i + dimensions];
        }
        partdist += alpha * alpha * rdiag[dimensions + j];
      }
      buffer.init_subtree(0, index, partdist, buffer.get_center_partsum(0, index, dimensions - 1));
    }
    debug_message_thread("Get %d new nodes\n", opts.initial_nodes_per_group);

    group.sync();

    clear_level<CG, block_size, levels, dimensions_per_level, max_nodes_per_level>(
        group, prefix_counter, group_shared_counter, buffer, 0, mu, rdiag, radius_squared_location,
        opts.tree_clear_opts, node_counter);
  }
}

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level, unsigned int start_point_dim>
void search(const std::array<std::array<float, levels * dimensions_per_level + start_point_dim>,
                             levels * dimensions_per_level + start_point_dim> &mu,
            const std::vector<std::array<enumi, start_point_dim>> &start_points,
            Opts<levels, dimensions_per_level, max_nodes_per_level> opts)
{
  typedef SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> SubtreeBuffer;

  assert(opts.start_point_dim == start_point_dim);

  constexpr unsigned int dimensions = dimensions_per_level * levels;
  constexpr unsigned int mu_n       = dimensions + start_point_dim;

  const unsigned int grid_size   = opts.thread_count / search_block_size;
  const unsigned int group_count = opts.thread_count / 32;
  const unsigned int group_size  = 32;
  assert(max_nodes_per_level >= opts.tree_clear_opts.max_subtree_paths * group_size);

  CudaPtr<unsigned char> buffer_mem =
      alloc(unsigned char, SubtreeBuffer::memory_size_in_bytes *group_count);
  CudaPtr<uint32_t> radius_mem             = alloc(uint32_t, 1);
  CudaPtr<enumf> device_mu                 = alloc(enumf, mu_n * mu_n);
  CudaPtr<enumf> device_rdiag              = alloc(enumf, mu_n);
  CudaPtr<unsigned long long> node_counter = alloc(unsigned long long, 5);
  CudaPtr<enumi> device_start_points       = alloc(enumi, start_points.size() * opts.start_point_dim);
  CudaPtr<unsigned int> processed_start_point_count = alloc(unsigned int, 1);

  const enumf radius                     = find_initial_radius(mu) * 1.01;
  const uint32_t radius_squared_location = float_to_int_order_preserving_bijection(radius * radius);
  std::unique_ptr<float[]> host_mu(new float[mu_n * mu_n]);
  std::unique_ptr<float[]> host_rdiag(new float[mu_n]);
  for (unsigned int i = 0; i < mu_n; ++i)
  {
    host_rdiag[i] = mu[i][i];
    for (unsigned int j = 0; j < mu_n; ++j)
    {
      host_mu[i * mu_n + j] = mu[i][j] / host_rdiag[i];
    }
    host_rdiag[i] = host_rdiag[i] * host_rdiag[i];
  }

  check(cudaMemcpy(device_mu.get(), host_mu.get(), mu_n * mu_n * sizeof(enumf),
                   cudaMemcpyHostToDevice));
  check(cudaMemcpy(device_rdiag.get(), host_rdiag.get(), mu_n * sizeof(enumf),
                   cudaMemcpyHostToDevice));
  check(cudaMemcpy(radius_mem.get(), &radius_squared_location, sizeof(uint32_t),
                   cudaMemcpyHostToDevice));
  check(cudaMemcpy(device_start_points.get(), start_points[0].data(),
                   start_point_dim * start_points.size() * sizeof(enumi), cudaMemcpyHostToDevice));

  std::cout << "started " << grid_size << " block with " << search_block_size << " threads each"
            << std::endl;

  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  search_kernel<<<dim3(grid_size), dim3(search_block_size)>>>(
      buffer_mem.get(), device_start_points.get(), processed_start_point_count.get(),
      start_points.size(), device_mu.get(), device_rdiag.get(), radius_mem.get(),
      node_counter.get(), opts);

  check(cudaDeviceSynchronize());
  check(cudaGetLastError());
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
            << std::endl;

  unsigned long long searched_nodes;
  uint32_t result_radius;
  check(cudaMemcpy(&searched_nodes, node_counter.get(), sizeof(unsigned long long),
                   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(&result_radius, radius_mem.get(), sizeof(uint32_t), cudaMemcpyDeviceToHost));
  std::cout << "searched nodes: " << searched_nodes << std::endl;
  std::cout << "result radius: " << sqrt(int_to_float_order_preserving_bijection(result_radius))
            << std::endl;
  print_performance_counter(&node_counter.get()[1]);
  print_performance_counter(&node_counter.get()[2]);
  print_performance_counter(&node_counter.get()[3]);
  print_performance_counter(&node_counter.get()[4]);
}
