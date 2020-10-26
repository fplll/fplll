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

#include "testdata.h"
#include "atomic.cuh"
#include "prefix.cuh"
#include "util.h"

typedef float enumf;

struct NodeCounter
{
  unsigned long long *counter;

  __device__ __host__ inline NodeCounter(unsigned long long *target) : counter(target) {}

  __device__ __host__ inline void count() { aggregated_atomic_inc(counter);
  }

  __device__ __host__ inline void perf_count(unsigned int c, unsigned long long time)
  {
    if (thread_id() == 0)
    atomic_add(&counter[c], time);
  }
};

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct SubtreeEnumerationBuffer;

__device__ __host__ inline enumf next_coeff(enumf coeff, const float center)
{
  const enumf rounded_center = round(center);
  coeff                      = 2 * rounded_center - coeff;
  if (center >= rounded_center)
  {
    return coeff + static_cast<int>(coeff >= rounded_center);
  }
  else
  {
    return coeff - static_cast<int>(coeff <= rounded_center);
  }
}

template <unsigned int maxdim> struct CudaEnumeration
{
  enumf x[maxdim];
  enumf partdist[maxdim];
  unsigned int center_partsum_begin[maxdim];
  // ! different to base enumeration of fplll, the second index is shifted !
  // _[i][j] contains inner product of i-th orthogonalized basis vector with B * (0, ..., 0, x[j +
  // 1], ... x[n])
  enumf center_partsums[maxdim][maxdim];
  enumf center[maxdim];
  const uint32_t *radius_squared_location;

  // row-major
  const enumf *mu;
  unsigned int ldmu;
  const enumf *rdiag;

  template <int kk, typename Callback>
  __device__ __host__ bool enumerate_recursive(Callback &, unsigned int &max_paths, NodeCounter& counter);

  template <int kk> __device__ __host__ bool is_enumeration_done() const;

  __device__ __host__ enumf get_radius_squared() const;

  template <unsigned int levels, unsigned int dimensions_per_level,
            unsigned int max_nodes_per_level>
  friend struct SubtreeEnumerationBuffer;
};

template <unsigned int maxdim>
__device__ __host__ inline enumf CudaEnumeration<maxdim>::get_radius_squared() const
{
  return int_to_float_order_preserving_bijection(*radius_squared_location);
}

template <unsigned int maxdim>
template <int kk>
__device__ __host__ inline bool CudaEnumeration<maxdim>::is_enumeration_done() const
{
  if constexpr (kk >= 0)
  {
    return isnan(x[kk]);
  }
  else
  {
    return false;
  }
}

/**
 * Searches the subtree of height kk + 1 using as root the values stored in this object. The
 * reference max_paths contains an integer that gives the maximal count of tree paths to search
 * (including tree paths that lead to nodes with too great partdist and are therefore cut). After
 * this is exceeded, the tree search is aborted but can be resumed later by calling
 * enumerate_recursive on this object. The function returns whether the subtree was completely
 * searched.
 */
template <unsigned int maxdim>
template <int kk, typename Callback>
__device__ __host__ __forceinline bool
CudaEnumeration<maxdim>::enumerate_recursive(Callback &callback, unsigned int &max_paths,
                                             NodeCounter& counter)
{
  static_assert(kk < static_cast<int>(maxdim),
                "Tree level count must be <= maximal enumeration dimension count");
  assert(max_paths >= 1);
  if constexpr (kk >= 0)
  {
    enumf alphak     = x[kk] - center[kk];
    enumf newdist = partdist[kk] + alphak * alphak * rdiag[kk];

    if (!(newdist <= get_radius_squared()))
    {
      x[kk] = NAN;
      return true;
    }

    if constexpr (kk == 0)
    {
      callback(x, newdist);
    }
    else
    {
      partdist[kk - 1] = newdist;

      for (int j = center_partsum_begin[kk]; j > kk - 1; --j)
      {
        center_partsums[kk - 1][j - 1] =
            center_partsums[kk - 1][j] - x[j] * mu[(kk - 1) * ldmu + j];
      }
      assert(!isnan(center_partsums[kk - 1][kk - 1]));

      if (center_partsum_begin[kk] > center_partsum_begin[kk - 1])
      {
        center_partsum_begin[kk - 1] = center_partsum_begin[kk];
      }

      center_partsum_begin[kk] = kk;
      center[kk - 1]           = center_partsums[kk - 1][kk - 1];
      if (isnan(x[kk - 1]))
      {
        x[kk - 1] = round(center[kk - 1]);
      }
    }

    while (true)
    {
      counter.count();
      bool is_done = enumerate_recursive<kk - 1, Callback>(callback, max_paths, counter);
      if (!is_done)
      {
        return false;
      }

      x[kk] = next_coeff(x[kk], center[kk]);

      enumf alphak2  = x[kk] - center[kk];
      enumf newdist2 = partdist[kk] + alphak2 * alphak2 * rdiag[kk];
      assert(!isnan(newdist2));

      if (max_paths == 1)
      {
        return false;
      }
      --max_paths;

      if (!(newdist2 <= get_radius_squared()))
      {
        x[kk] = NAN;
        return true;
      }

      if constexpr (kk == 0)
      {
        callback(x, newdist2);
      }
      else
      {
        partdist[kk - 1] = newdist2;

        center_partsums[kk - 1][kk - 1] =
            center_partsums[kk - 1][kk - 1 + 1] - x[kk - 1 + 1] * mu[(kk - 1) * ldmu + kk - 1 + 1];
        assert(!isnan(center_partsums[kk - 1][kk - 1]));

        if (kk > center_partsum_begin[kk - 1])
          center_partsum_begin[kk - 1] = kk;

        center[kk - 1] = center_partsums[kk - 1][kk - 1];
        x[kk - 1]      = round(center[kk - 1]);
      }
    }
  }
  return true;
}

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct SubtreeEnumerationBuffer
{
private:
  // shape [levels, dimensions_per_level, max_nodes_per_level]
  enumf *enumeration_x;
  // shape [levels, dimensions_per_level, max_nodes_per_level]
  enumf *coefficients;
  // shape [levels, dimensions, max_nodes_per_level], of subtree root
  enumf *center_partsum;
  // shape [levels, max_nodes_per_level], of subtree root
  enumf *partdist;
  // shape [levels, max_nodes_per_level]
  unsigned int *parent_indices;

  // shape [levels]
  unsigned int *open_node_count;

  constexpr static unsigned int dimensions = levels * dimensions_per_level;

  constexpr static unsigned int enumeration_x_memory_size_in_bytes =
      sizeof(enumf) * levels * dimensions_per_level * max_nodes_per_level;

  constexpr static unsigned int coefficient_memory_size_in_bytes =
      sizeof(enumf) * levels * dimensions_per_level * max_nodes_per_level;

  constexpr static unsigned int center_partsum_size_in_bytes =
      sizeof(enumf) * levels * dimensions * max_nodes_per_level;

  constexpr static unsigned int partdist_memory_size_in_bytes =
      sizeof(enumf) * levels * max_nodes_per_level;

  constexpr static unsigned int parent_indices_memory_size_in_bytes =
      sizeof(unsigned int) * levels * max_nodes_per_level;

  constexpr static unsigned int open_node_count_memory_size_in_bytes =
      sizeof(unsigned int) * levels;

public:
  __device__ __host__ inline SubtreeEnumerationBuffer(unsigned char *memory)
      : enumeration_x(reinterpret_cast<enumf *>(memory)),
        coefficients(reinterpret_cast<enumf *>(memory + enumeration_x_memory_size_in_bytes)),
        center_partsum(reinterpret_cast<enumf *>(memory + enumeration_x_memory_size_in_bytes +
                                                 coefficient_memory_size_in_bytes)),
        partdist(reinterpret_cast<enumf *>(memory + enumeration_x_memory_size_in_bytes +
                                           coefficient_memory_size_in_bytes +
                                           center_partsum_size_in_bytes)),
        parent_indices(reinterpret_cast<unsigned int *>(
            memory + enumeration_x_memory_size_in_bytes + coefficient_memory_size_in_bytes +
            center_partsum_size_in_bytes + partdist_memory_size_in_bytes)),
        open_node_count(reinterpret_cast<unsigned int *>(
            memory + enumeration_x_memory_size_in_bytes + coefficient_memory_size_in_bytes +
            center_partsum_size_in_bytes + partdist_memory_size_in_bytes +
            parent_indices_memory_size_in_bytes))
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
      enumeration_x_memory_size_in_bytes + coefficient_memory_size_in_bytes +
      center_partsum_size_in_bytes + partdist_memory_size_in_bytes +
      parent_indices_memory_size_in_bytes + open_node_count_memory_size_in_bytes;

  static_assert(memory_size_in_bytes < std::numeric_limits<unsigned int>::max(),
                "Requires more memory than indexable with unsigned int");

  __device__ __host__ CudaEnumeration<dimensions_per_level>
  get_enumeration(unsigned int tree_level, unsigned int index, const float *mu, const unsigned int ldmu, const float *rdiag,
                  const uint32_t *radius_squared_location)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    CudaEnumeration<dimensions_per_level> result;

    const unsigned int offset_kk = (levels - tree_level - 1) * dimensions_per_level;

    result.mu                      = &mu[offset_kk * ldmu + offset_kk];
    result.ldmu                    = dimensions;
    result.rdiag                   = &rdiag[offset_kk];
    result.radius_squared_location = radius_squared_location;

    for (unsigned int i = 0; i < dimensions_per_level; ++i)
    {
      result.x[i] = enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                                  i * max_nodes_per_level + index];

      result.center_partsum_begin[i] = dimensions_per_level - 1;

      const enumf center_partsum_i =
          center_partsum[tree_level * dimensions * max_nodes_per_level +
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

  __device__ __host__ void set_enumeration(unsigned int tree_level, unsigned int index,
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

  __device__ __host__ void init_subtree(unsigned int tree_level, unsigned int index,
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

  __device__ __host__ void set_center_partsum(unsigned int tree_level, unsigned int index,
                                              unsigned int orth_basis_index, enumf value)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    assert(orth_basis_index < dimensions);
    center_partsum[tree_level * dimensions * max_nodes_per_level +
                   orth_basis_index * max_nodes_per_level + index] = value;
  }

  __device__ __host__ enumf get_center_partsum(unsigned int tree_level, unsigned int index,
                                               unsigned int orth_basis_index)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    assert(orth_basis_index < dimensions);
    return center_partsum[tree_level * dimensions * max_nodes_per_level +
                          orth_basis_index * max_nodes_per_level + index];
  }

  __device__ __host__ enumf get_coefficient(unsigned int tree_level, unsigned int index,
                                            unsigned int coordinate)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    assert(coordinate < dimensions_per_level);
    return coefficients[tree_level * dimensions_per_level * max_nodes_per_level +
                        coordinate * max_nodes_per_level + index];
  }

  __device__ __host__ void set_coefficient(unsigned int tree_level, unsigned int index,
                                           unsigned int coordinate, enumf value)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    assert(coordinate < dimensions_per_level);
    coefficients[tree_level * dimensions_per_level * max_nodes_per_level +
                 coordinate * max_nodes_per_level + index] = value;
  }

  __device__ __host__ unsigned int get_parent_index(unsigned int tree_level, unsigned int index)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    return parent_indices[tree_level * max_nodes_per_level + index];
  }

  __device__ __host__ inline unsigned int get_node_count(unsigned int tree_level)
  {
    assert(tree_level < levels);
    return open_node_count[tree_level];
  }

  __device__ __host__ inline unsigned int add_subtree(unsigned int tree_level,
                                                      unsigned int parent_node_index)
  {
    assert(tree_level < levels);
    const unsigned int new_task_index = aggregated_atomic_inc(&open_node_count[tree_level]);
    assert(new_task_index < max_nodes_per_level);
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

    enumf coefficients_tmp[dimensions_per_level];
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
  const enumf *mu;
  unsigned int ldmu;
  uint32_t *radius_squared_location;
  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer;

  __device__ __host__ void operator()(const float *x, float squared_norm);
};

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ __forceinline void
ProcessLeafCallback<levels, dimensions_per_level, max_nodes_per_level>::operator()(
    const float *x, float squared_norm)
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
    if (TRACE)
    {
      printf("Squared norm %f: ", squared_norm);
      float coefficient;
      for (unsigned int i = 0; i < dimensions_per_level; ++i)
      {
        coefficient = x[i];
        printf("%f, ", coefficient);
      }
      unsigned int index = parent_index;
      for (int j = levels - 1; j > 0; --j)
      {
        for (unsigned int i = 0; i < dimensions_per_level; ++i)
        {
          coefficient = buffer.get_coefficient(j, index, i);
          printf("%f, ", coefficient);
        }
        index = buffer.get_parent_index(j, index);
      }
      printf("\n");
    }
  }
}

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct AddToTreeCallback
{
  unsigned int level;
  unsigned int parent_index;
  const enumf *mu;
  unsigned int ldmu;
  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer;
  NodeCounter &counter;

  __device__ __host__ void operator()(const float *x, float squared_norm);
};

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ __forceinline void
AddToTreeCallback<levels, dimensions_per_level, max_nodes_per_level>::operator()(const float *x,
                                                                                 float squared_norm)
{
  assert(level > 0);

  unsigned long long begin     = time();

  const unsigned int new_index = buffer.add_subtree(level, parent_index);
  unsigned int kk_offset       = (levels - level - 1) * dimensions_per_level;
  for (unsigned int i = 0; i < kk_offset + dimensions_per_level; ++i)
  {
    enumf center_partsum = buffer.get_center_partsum(level - 1, parent_index, i);
    for (unsigned int j = 0; j < dimensions_per_level; ++j)
    {
      center_partsum -= x[j] * mu[i * ldmu + j + dimensions_per_level + kk_offset];
    }
    assert(!isnan(center_partsum));
    buffer.set_center_partsum(level, new_index, i, center_partsum);
  }
  unsigned int subtree_root_expand_vector_index = kk_offset + dimensions_per_level;
  enumf center = buffer.get_center_partsum(level, new_index, kk_offset + dimensions_per_level - 1);
  for (unsigned int j = 0; j < dimensions_per_level; ++j)
  {
    buffer.set_coefficient(level, new_index, j, x[j]);
  }
  buffer.init_subtree(level, new_index, squared_norm, center);

  unsigned long long end = time();
  counter.perf_count(1, end - begin);
}

template <typename CG, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ void
do_search_step(CG &group,
               SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
               unsigned int level, const float *mu, unsigned int ldmu, const float *rdiag,
               uint32_t *radius_squared_location, unsigned int max_subtree_paths, NodeCounter& counter)
{
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  unsigned int max_paths = max_subtree_paths;
  if (active)
  {
    CudaEnumeration<dimensions_per_level> enumeration =
        buffer.get_enumeration(level, index, mu, ldmu, rdiag, radius_squared_location);

    unsigned long long begin = time();
    if (level == levels - 1)
    {
      typedef ProcessLeafCallback<levels, dimensions_per_level, max_nodes_per_level> CallbackT;
      CallbackT callback = {level + 1, index, mu, ldmu, radius_squared_location, buffer};
      enumeration.template enumerate_recursive<dimensions_per_level - 1, CallbackT>(callback, max_paths, counter);
    }
    else
    {
      typedef AddToTreeCallback<levels, dimensions_per_level, max_nodes_per_level> CallbackType;
      CallbackType callback = {level + 1, index, mu, ldmu, buffer, counter};
      enumeration.template enumerate_recursive<dimensions_per_level - 1, CallbackType>(callback, max_paths, counter);
    }
    unsigned long long end = time();
    counter.perf_count(2, end - begin);

    buffer.set_enumeration(level, index, enumeration);
  }
}

template <typename CG, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ inline void get_done_subtree_count(
    CG &group, unsigned int *shared_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
    unsigned int level, const float *mu, unsigned int ldmu, const float *rdiag, const uint32_t *radius_square_location)
{
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  if (active)
  {
    bool is_done = buffer.get_enumeration(level, index, mu, ldmu, rdiag, radius_square_location)
                       .template is_enumeration_done<dimensions_per_level - 1>();
    if (is_done)
    {
      aggregated_atomic_inc(shared_counter);
    }
  }
}

template <typename CG, unsigned int block_size, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void
do_cleanup_step(CG &group, PrefixCounter<CG, block_size> &prefix_counter,
                SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
                unsigned int level, const float *mu, unsigned int ldmu, const float *rdiag,
                const uint32_t *radius_square_location)
{
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  bool is_done = buffer.get_enumeration(level, index, mu, ldmu, rdiag, radius_square_location)
                     .template is_enumeration_done<dimensions_per_level - 1>();
  buffer.filter_nodes(group, prefix_counter, level, index, !is_done, active_thread_count);
}

template <typename CG, unsigned int block_size, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void
clear_level(CG &group, PrefixCounter<CG, block_size> &prefix_counter, unsigned int *shared_counter,
            SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
            int level, const float *mu, unsigned int ldmu, const float *rdiag, uint32_t *radius_square_location,
            unsigned int max_subtree_paths, NodeCounter counter)
{
  unsigned long long begin_outer = time();
  while (level >= 0)
  {
    assert(all_threads_eq(group, level, shared_counter));
    if (level + 1 < levels)
    {
      if (buffer.get_node_count(level) > 0)
      {
        const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
        // create as many children as fit into the children buffer (and only as long as we can generate using enough parallel threads)
        while (true)
        {
          do_search_step(group, buffer, level, mu, ldmu, rdiag, radius_square_location, max_subtree_paths,
                         counter);
          if (group.thread_rank() == 0)
          {
            *shared_counter = 0;
          }
          group.sync();
          get_done_subtree_count(group, shared_counter, buffer, level, mu, ldmu, rdiag,
                                 radius_square_location);
          group.sync();
          
          debug_message_thread(
                "Worked on level %d, next level points are %d, %d nodes of current working pool "
                   "(%d) are done\n",
                   level, buffer.get_node_count(level + 1), *shared_counter, active_thread_count);

          if (buffer.get_node_count(level + 1) >= max_nodes_per_level / 2)
          {
            break;
          }
          else if (*shared_counter >= active_thread_count / 2)
          {
            break;
          }
          group.sync();
        }
        ++level;
      }
      else
      {
        --level;
        if (level >= 0)
        {
          group.sync();
          do_cleanup_step(group, prefix_counter, buffer, level, mu, ldmu, rdiag, radius_square_location);
          group.sync();

          debug_message_thread("Cleaned up level %d, has now %d nodes\n",
                   level, buffer.get_node_count(level));

          group.sync();
        }
      }
    }
    else
    {
      while (buffer.get_node_count(level) > 0)
      {
        do_search_step(group, buffer, level, mu, ldmu, rdiag, radius_square_location, max_subtree_paths,
                       counter);
        do_cleanup_step(group, prefix_counter, buffer, level, mu, ldmu, rdiag, radius_square_location);
      }
      --level;
      if (level >= 0)
      {
        group.sync();
        do_cleanup_step(group, prefix_counter, buffer, level, mu, ldmu, rdiag, radius_square_location);
        group.sync();
      }
    }
  }
  unsigned long long end_outer = time();
  counter.perf_count(3, end_outer - begin_outer);
}

template <unsigned int block_size, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__global__ void __launch_bounds__(512, 1)
    search_kernel(unsigned char *buffer_memory, const enumf *start_points,
                  unsigned int *processed_start_point_counter, unsigned int start_point_dim, unsigned int start_point_count, const enumf *mu, const enumf *rdiag,
                  uint32_t *radius_squared_location, unsigned int max_subtree_paths,
                  unsigned long long *counter)
{
  typedef cooperative_groups::thread_block_tile<32> CG;
  typedef SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> SubtreeBuffer;

  constexpr unsigned int dimensions = dimensions_per_level * levels;

  const unsigned int group_count_per_block = block_size / 32;
  extern __shared__ unsigned char
      shared_mem[PrefixCounter<CG, block_size>::shared_mem_size_in_bytes() +
                 group_count_per_block * sizeof(unsigned int)];

  CG group = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());
  const unsigned int group_id = thread_id() / 32;
  const unsigned int group_id_in_block = thread_id_in_block() / 32;
  PrefixCounter<CG, block_size> prefix_counter /*(shared_mem + group_count_per_block * sizeof(unsigned int))*/;
  unsigned int *shared_counter = reinterpret_cast<unsigned int *>(shared_mem + group_id_in_block * sizeof(unsigned int));

  NodeCounter node_counter(counter);
  constexpr unsigned int nodes_per_group = 4;
  assert(nodes_per_group <= group.size());

  const unsigned int ldmu = dimensions + start_point_dim;

  SubtreeBuffer buffer(buffer_memory + group_id * SubtreeBuffer::memory_size_in_bytes);

  while (true)
  {
    group.sync();
    if (group.thread_rank() == 0)
    {
      *shared_counter = atomic_add(processed_start_point_counter, nodes_per_group);
    }
    buffer.init(group);
    group.sync();

    if (*shared_counter >= start_point_count)
    {
      return;
    }
    const unsigned int start_point_index = *shared_counter + group.thread_rank();
    const bool active =
        group.thread_rank() < nodes_per_group && start_point_index < start_point_count;

    if (active)
    {
      const unsigned int index = buffer.add_subtree(0, start_point_index);
      for (unsigned int i = 0; i < dimensions; ++i)
      {
        enumf partsum = 0;
        for (unsigned int j = 0; j < start_point_dim; ++j)
        {
          partsum -=
              start_points[start_point_index * start_point_dim + j] * mu[i * ldmu + dimensions + j];
        }
        buffer.set_center_partsum(0, index, i, partsum);
      }
      enumf partdist = 0;
      for (unsigned int j = 0; j < start_point_dim; ++j)
      {
        enumf coefficient = start_points[start_point_index * start_point_dim + j];
        partdist += coefficient * coefficient * rdiag[dimensions + j];
      }
      buffer.init_subtree(0, index, partdist, buffer.get_center_partsum(0, index, dimensions - 1));
    }

    group.sync();

    clear_level<CG, block_size, levels, dimensions_per_level, max_nodes_per_level>(
        group, prefix_counter, shared_counter, buffer, 0, mu, ldmu, rdiag, radius_squared_location,
        max_subtree_paths, node_counter);
  }
}

template <unsigned int block_size, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level, unsigned int start_point_dim>
void search(const std::array<std::array<enumf, levels * dimensions_per_level + start_point_dim>,
                             levels * dimensions_per_level + start_point_dim> &mu,
            const std::vector < std::array < enumf, start_point_dim>>& start_points,
            unsigned int max_subtree_paths, unsigned int grid_size)
{
  typedef SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> SubtreeBuffer;

  constexpr unsigned int dimensions = dimensions_per_level * levels;
  constexpr unsigned int mu_n       = dimensions + start_point_dim;

  const unsigned int group_count = grid_size * block_size / 32;
  const unsigned int group_size  = 32;
  assert(max_nodes_per_level >= max_subtree_paths * group_size);

  CudaPtr<unsigned char> buffer_mem =
      alloc(unsigned char, SubtreeBuffer::memory_size_in_bytes *group_count);
  CudaPtr<uint32_t> radius_mem        = alloc(uint32_t, 1);
  CudaPtr<enumf> device_mu            = alloc(enumf, mu_n * mu_n);
  CudaPtr<enumf> device_rdiag         = alloc(enumf, mu_n);
  CudaPtr<unsigned long long> counter = alloc(unsigned long long, 4);
  CudaPtr<enumf> device_start_points  = alloc(enumf, start_points.size() * start_point_dim);
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
                   start_point_dim * start_points.size() * sizeof(enumf), cudaMemcpyHostToDevice));

  std::cout << "started " << grid_size << " block with " << block_size << " threads each"
            << std::endl;

  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  search_kernel<block_size, levels, dimensions_per_level, max_nodes_per_level>
      <<<dim3(grid_size), dim3(block_size)>>>(
          buffer_mem.get(), device_start_points.get(), processed_start_point_count.get(),
          start_point_dim, start_points.size(), device_mu.get(), device_rdiag.get(), radius_mem.get(),
          max_subtree_paths, counter.get());

  check(cudaDeviceSynchronize());
  check(cudaGetLastError());
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
            << std::endl;

  unsigned long long searched_nodes;
  uint32_t result_radius;
  check(cudaMemcpy(&searched_nodes, counter.get(), sizeof(unsigned long long),
                   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(&result_radius, radius_mem.get(), sizeof(uint32_t),
                   cudaMemcpyDeviceToHost));
  std::cout << "searched nodes: " << searched_nodes << std::endl;
  std::cout << "result radius: " << sqrt(int_to_float_order_preserving_bijection(result_radius))
            << std::endl;
  print_performance_counter(&counter.get()[1]);
  print_performance_counter(&counter.get()[2]);
  print_performance_counter(&counter.get()[3]);
}
