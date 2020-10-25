#pragma once
#include "cooperative_groups.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <memory>
#include <assert.h>

#include "util.h"
#include "atomic.cuh"
#include "prefix.cuh"
#include "testdata.h"

namespace hybrid_enum
{

typedef float enumf;

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

template<unsigned int maxdim>
struct CudaEnumeration
{
  enumf x[maxdim];
  enumf partdist[maxdim];
  unsigned int center_partsum_begin[maxdim];
  // ! different to base enumeration of fplll, the second index is shifted !
  // _[i][j] contains inner product of i-th orthogonalized basis vector with B * (0, ..., 0, x[j + 1], ... x[n])
  enumf center_partsums[maxdim][maxdim];
  enumf center[maxdim];
  enumf radius_squared;

  // row-major
  const enumf *mu;
  unsigned int ldmu;
  const enumf *rdiag;

  template <int kk, typename Callback> __device__ __host__ bool enumerate_recursive(Callback &, unsigned int& max_paths);

  template <int kk> __device__ __host__ bool is_enumeration_done() const;

  template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
  friend struct SubtreeEnumerationBuffer;
};

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
 * Searches the subtree of height kk + 1 using as root the values stored in this object. The reference max_paths
 * contains an integer that gives the maximal count of tree paths to search (including tree paths that lead to nodes
 * with too great partdist and are therefore cut). After this is exceeded, the tree search is aborted but can be
 * resumed later by calling enumerate_recursive on this object. The function returns whether the subtree was completely
 * searched.
 */
template <unsigned int maxdim>
template <int kk, typename Callback>
__device__ __host__ inline bool
CudaEnumeration<maxdim>::enumerate_recursive(Callback &callback, unsigned int& max_paths)
{
  static_assert(kk < static_cast<int>(maxdim), "Tree level count must be <= maximal enumeration dimension count");
  assert(max_paths >= 1);
  if constexpr (kk >= 0)
  {
    enumf alphak  = x[kk] - center[kk];
    enumf newdist = partdist[kk] + alphak * alphak * rdiag[kk];

    if (newdist > radius_squared)
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
        center_partsums[kk - 1][j - 1] =
            center_partsums[kk - 1][j] - x[j] * mu[(kk - 1) * ldmu + j];

      if (center_partsum_begin[kk] > center_partsum_begin[kk - 1])
        center_partsum_begin[kk - 1] = center_partsum_begin[kk];

      center_partsum_begin[kk] = kk;
      center[kk - 1]           = center_partsums[kk - 1][kk - 1];
      if (isnan(x[kk - 1]))
      {
        x[kk - 1] = round(center[kk - 1]);
      }
    }

    while (true)
    {
      bool is_done = enumerate_recursive<kk - 1, Callback>(callback, max_paths);
      if (!is_done)
      {
        return false;
      }

      x[kk] = next_coeff(x[kk], center[kk]);

      enumf alphak2  = x[kk] - center[kk];
      enumf newdist2 = partdist[kk] + alphak2 * alphak2 * rdiag[kk];

      if (max_paths == 1)
      {
        return false;
      }
      --max_paths;

      if (newdist2 > radius_squared)
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

        center_partsums[kk - 1][kk - 1] = center_partsums[kk - 1][kk - 1 + 1] - x[kk - 1 + 1] * mu[(kk - 1) * ldmu + kk - 1 + 1];

        if (kk > center_partsum_begin[kk - 1])
          center_partsum_begin[kk - 1] = kk;

        center[kk - 1] = center_partsums[kk - 1][kk - 1];
        x[kk - 1] = round(center[kk - 1]);
      }
    }
  }
  else
  {
    return true;
  }
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
      center_partsum_size_in_bytes +
      partdist_memory_size_in_bytes + parent_indices_memory_size_in_bytes +
      open_node_count_memory_size_in_bytes;

  static_assert(memory_size_in_bytes < std::numeric_limits<unsigned int>::max(),
                "Requires more memory than indexable with unsigned int");

  __device__ __host__ CudaEnumeration<dimensions_per_level>
  get_enumeration(unsigned int tree_level, unsigned int index, const float *mu, const float *rdiag,
                  const float radius_squared)
  {
    assert(tree_level < levels);
    assert(index < max_nodes_per_level);
    CudaEnumeration<dimensions_per_level> result;

    const unsigned int offset_kk = (levels - tree_level - 1) * dimensions_per_level;

    result.mu                    = mu + offset_kk * dimensions + offset_kk;
    result.ldmu                  = dimensions;
    result.rdiag                 = rdiag + offset_kk;
    result.radius_squared        = radius_squared;

    for (unsigned int i = 0; i < dimensions_per_level; ++i)
    {
      result.x[i] = enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
                                 i * max_nodes_per_level + index];
      result.center_partsum_begin[i] = dimensions_per_level - 1;
      result.center_partsums[i][dimensions_per_level - 1] =
          center_partsum[tree_level * dimensions * max_nodes_per_level +
                         (offset_kk + i) * max_nodes_per_level +
                         index];
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

  __device__ __host__ void init_subtree(unsigned int tree_level, unsigned int index, enumf parent_partdist, enumf center) {
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
    return center_partsum[tree_level * dimensions * max_nodes_per_level + orth_basis_index * max_nodes_per_level + index];
  }

  __device__ __host__ enumf get_coefficient(unsigned int tree_level, unsigned int index, unsigned int coordinate)
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

  __device__ __host__ inline unsigned int
  add_subtree(unsigned int tree_level, unsigned int parent_node_index)
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
                            unsigned int tree_level, unsigned int old_index,
                            bool keep_this_thread_task, unsigned int active_thread_count)
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
      partdist_tmp = partdist[tree_level * max_nodes_per_level + old_index];
      parent_index_tmp = parent_indices[tree_level * max_nodes_per_level + old_index];
      for (unsigned int i = 0; i < dimensions_per_level; ++i)
      {
        coefficients_tmp[i] = enumeration_x[tree_level * dimensions_per_level * max_nodes_per_level +
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
      partdist[tree_level * max_nodes_per_level + new_index] = partdist_tmp;
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
  constexpr static unsigned int ldmu = levels * dimensions_per_level;

  unsigned int level;
  unsigned int parent_index;
  const enumf *mu;
  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer;

  __device__ __host__ void operator()(const float *x, float squared_norm);
};

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void ProcessLeafCallback<levels, dimensions_per_level, max_nodes_per_level>::operator()(const float *x,
                                                                  float squared_norm)
{
  printf("Squared norm %f: ", squared_norm);
  float coefficient;
  for (unsigned int i = 0; i < dimensions_per_level; ++i)
  {
    coefficient = x[i];
    printf("%f, ", coefficient);
  }
  unsigned int index = parent_index;
  for (int j = levels - 1; j > 0 ; --j)
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

template <unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
struct AddToTreeCallback
{
  constexpr static unsigned int ldmu = levels * dimensions_per_level;

  unsigned int level;
  unsigned int parent_index;
  const enumf *mu;
  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer;

  __device__ __host__ void operator()(const float *x, float squared_norm);
};

template <
    unsigned int levels, unsigned int dimensions_per_level,
    unsigned int max_nodes_per_level>
__device__ __host__ inline void
AddToTreeCallback<levels, dimensions_per_level, max_nodes_per_level>::operator()(const float *x, float squared_norm)
{
  assert(level > 0);
  const unsigned int new_index = buffer.add_subtree(level, parent_index);
  unsigned int kk_offset       = (levels - level - 1) * dimensions_per_level;
  for (unsigned int i = 0; i < kk_offset + dimensions_per_level; ++i)
  {
    enumf center_partsum = buffer.get_center_partsum(level - 1, parent_index, i);
    for (unsigned int j = 0; j < dimensions_per_level; ++j)
    {
      center_partsum -= x[j] * mu[i * ldmu + j + dimensions_per_level + kk_offset];
    }
    buffer.set_center_partsum(level, new_index, i, center_partsum);
  }
  unsigned int subtree_root_expand_vector_index = kk_offset + dimensions_per_level;
  enumf center = buffer.get_center_partsum(level, new_index, kk_offset + dimensions_per_level - 1);
  for (unsigned int j = 0; j < dimensions_per_level; ++j)
  {
    buffer.set_coefficient(level, new_index, j, x[j]);
  }
  buffer.init_subtree(level, new_index, squared_norm, center);
}

template <typename CG, unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ void
do_search_step(CG &group,
               SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
               unsigned int level, const float *mu, const float *rdiag, float radius_square, unsigned int max_subtree_paths)
{
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  unsigned int max_paths = max_subtree_paths;
  if (active)
  {
    CudaEnumeration<dimensions_per_level> enumeration =
        buffer.get_enumeration(level, index, mu, rdiag, radius_square);
    if (level == levels - 1)
    {
      typedef ProcessLeafCallback<levels, dimensions_per_level, max_nodes_per_level> CallbackT;
      CallbackT callback = {level + 1, index, mu, buffer};
      enumeration.template enumerate_recursive<dimensions_per_level - 1, CallbackT>(callback,
                                                                                    max_paths);
    }
    else
    {
      typedef AddToTreeCallback<levels, dimensions_per_level, max_nodes_per_level> CallbackType;
      CallbackType callback = {level + 1, index, mu, buffer};
      enumeration.template enumerate_recursive<dimensions_per_level - 1, CallbackType>(callback,
                                                                                       max_paths);
    }
    buffer.set_enumeration(level, index, enumeration);
  }
}

template <typename CG, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void get_done_subtree_count(CG& group,
    unsigned int *shared_counter,
    SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
    unsigned int level, const float* mu, const float* rdiag, const float radius_square)
{
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  bool is_done = buffer.get_enumeration(level, index, mu, rdiag, radius_square)
                     .template is_enumeration_done<dimensions_per_level - 1>();
  if (is_done)
  {
    aggregated_atomic_inc(shared_counter);
  }
}

template <typename CG, unsigned int block_size, unsigned int levels, unsigned int dimensions_per_level,
          unsigned int max_nodes_per_level>
__device__ __host__ inline void
do_cleanup_step(CG &group, PrefixCounter<CG, block_size>& prefix_counter,
                SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
                unsigned int level, const float *mu, const float *rdiag, float radius_square)
{
  const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
  const unsigned int index =
      buffer.get_node_count(level) - active_thread_count + group.thread_rank();
  const bool active = index < buffer.get_node_count(level);

  bool is_done = buffer.get_enumeration(level, index, mu, rdiag, radius_square)
                     .template is_enumeration_done<dimensions_per_level - 1>();
  buffer.filter_nodes(group, prefix_counter, level, index, !is_done, active_thread_count);
}

template <typename CG, unsigned int block_size, unsigned int levels,
          unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
__device__ __host__ inline void clear_level(CG &group, PrefixCounter<CG, block_size> &prefix_counter,
                unsigned int* shared_counter,
                SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> &buffer,
            unsigned int level, const float *mu, const float *rdiag, float radius_square,
            unsigned int max_subtree_paths)
{
  if (level < levels - 1)
  {
    while (buffer.get_node_count(level) > 0)
    {
      const unsigned int active_thread_count = min(buffer.get_node_count(level), group.size());
      while (true)
      {
        do_search_step(group, buffer, level, mu, rdiag, radius_square, max_subtree_paths);
        *shared_counter = 0;
        group.sync();
        get_done_subtree_count(group, shared_counter, buffer, level, mu, rdiag, radius_square);
        group.sync();
        if (buffer.get_node_count(level + 1) + 10 >= max_nodes_per_level)
        {
          break;
        }
        else if (*shared_counter >= active_thread_count / 2)
        {
          break;
        }
        group.sync();
      }
      clear_level(group, prefix_counter, shared_counter, buffer, level + 1, mu, rdiag,
                  radius_square, max_subtree_paths);
      do_cleanup_step(group, prefix_counter, buffer, level, mu, rdiag, radius_square);
    }
  }
  else
  {
    while (buffer.get_node_count(level) > 0)
    {
      do_search_step(group, buffer, level, mu, rdiag, radius_square, max_subtree_paths);
      do_cleanup_step(group, prefix_counter, buffer, level, mu, rdiag, radius_square);
    }
  }
}
    
inline void cpu_test()
{
  constexpr unsigned int levels               = 4;
  constexpr unsigned int dimensions_per_level = 5;
  constexpr unsigned int dimensions           = levels * dimensions_per_level;
  constexpr unsigned int max_nodes_per_level  = 200000;

  const std::array<std::array<float, dimensions>, dimensions>& host_mu = test_mu_small;

  single_thread group;
  PrefixCounter<single_thread, 1> prefix_counter;

  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> buffer(
      new unsigned char[SubtreeEnumerationBuffer<levels, dimensions_per_level,
                                                 max_nodes_per_level>::memory_size_in_bytes]);
  buffer.init(group);

  std::unique_ptr<float[]> mu(new float[dimensions * dimensions]);
  std::unique_ptr<float[]> rdiag(new float[dimensions]);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    rdiag[i] = host_mu[i][i];
    for (unsigned int j = 0; j < dimensions; ++j)
    {
      mu[i * dimensions + j] = host_mu[i][j] / rdiag[i];
    }
    rdiag[i] = rdiag[i] * rdiag[i];
  }
  float radius             = find_initial_radius(host_mu) * 1.1;

  const unsigned int index = buffer.add_subtree(0, 0);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    buffer.set_center_partsum(0, index, i, 0);
  }
  buffer.init_subtree(0, index, 0, 0);

  unsigned int counter = 0;
  clear_level<single_thread, 1, levels, dimensions_per_level, max_nodes_per_level>(group, prefix_counter, &counter, buffer, 0, mu.get(), rdiag.get(), radius * radius, 5);
}

inline void cpu_test4d()
{
  constexpr unsigned int levels               = 1;
  constexpr unsigned int dimensions_per_level = 4;
  constexpr unsigned int dimensions           = levels * dimensions_per_level;
  constexpr unsigned int max_nodes_per_level  = 100;

  std::array<std::array<float, dimensions>, dimensions> host_mu = {
      {{3, 2, 4, 3}, {0, 3, 4, 2}, {0, 0, 3, 3}, {0, 0, 0, 2}}};

  single_thread group;
  PrefixCounter<single_thread, 1> prefix_counter;

  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> buffer(
      new unsigned char[SubtreeEnumerationBuffer<levels, dimensions_per_level,
                                                 max_nodes_per_level>::memory_size_in_bytes]);
  buffer.init(group);

  std::unique_ptr<float[]> mu(new float[dimensions * dimensions]);
  std::unique_ptr<float[]> rdiag(new float[dimensions]);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    rdiag[i] = host_mu[i][i];
    for (unsigned int j = 0; j < dimensions; ++j)
    {
      mu[i * dimensions + j] = host_mu[i][j] / rdiag[i];
    }
    rdiag[i] = rdiag[i] * rdiag[i];
  }
  float radius = 3.01;

  const unsigned int index = buffer.add_subtree(0, 0);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    buffer.set_center_partsum(0, index, i, 0);
  }
  buffer.init_subtree(0, index, 0, 0);

  unsigned int counter = 0;
  clear_level<single_thread, 1, levels, dimensions_per_level, max_nodes_per_level>(
      group, prefix_counter, &counter, buffer, 0, mu.get(), rdiag.get(), radius * radius, 5);
}

inline void test() { cpu_test(); }
}
