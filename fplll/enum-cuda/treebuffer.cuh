#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <assert.h>

#include "atomic.cuh"
#include "prefix.cuh"
#include "util.h"

/**
 * Stores a set of enumeration tree nodes whose subtrees have to be processed, ordered by tree level
 * in a way that maximizes memory coalescing.
 *
 * For general information about how the tree search works, refer to the readme.
 *
 * Template Parameters:
 * levels               Count of tree levels, excluding the root level
 *
 * dimension            Count of coefficients per vector; For standard enumeration, this is equal to
 *                      levels, but to distribute work among blocks, we might not start on the root level of the
 *                      enumeration tree
 *
 * max_nodes_per_level  Maximal count of open nodes on a tree level; Should be >= 2 * group size
 */
template <typename FL, unsigned int levels, unsigned int max_nodes_per_level,
          unsigned int dimension>
struct OpenExpansionTreeNodeBuffer
{
private:
  // shape [levels, dimension, max_nodes_per_level]; starting points are at level 0, tree leafs are
  // at level (levels - 1)
  FL *__restrict__ coefficients;
  // shape [levels, max_nodes_per_level]
  CenterToOutIterator *__restrict__ children_iterator;

  // shape [levels, max_nodes_per_level], contains the squared norm of the point projected into the
  // subspace perpendicular to the first k lattice vectors
  FL *__restrict__ partdist;
  // shape [levels, levels, max_nodes_per_level]; [k, j, i] contains the scalar product of the i-th
  // point at level k with the j-th orthogonalized basis vector
  FL *__restrict__ center_partsums;

  // shape [levels]
  unsigned int *__restrict__ open_node_count;

  constexpr static unsigned int coefficient_memory_size_in_bytes =
      sizeof(FL) * dimension * max_nodes_per_level * levels;

  constexpr static unsigned int children_iterator_size_in_bytes =
      sizeof(CenterToOutIterator) * max_nodes_per_level * levels;

  constexpr static unsigned int partdist_memory_size_in_bytes =
      sizeof(FL) * levels * max_nodes_per_level;

  constexpr static unsigned int center_partsums_memory_size_in_bytes =
      sizeof(FL) * levels * levels * max_nodes_per_level;

  constexpr static unsigned int open_node_count_memory_size_in_bytes =
      sizeof(unsigned int) * levels;

public:
  static_assert(
      dimension >= levels,
      "Dimension of coefficient vectors must be at least as great as the count of tree levels");

  __device__ __host__ inline OpenExpansionTreeNodeBuffer(unsigned char *memory)
      : coefficients(reinterpret_cast<FL *>(memory)),
        children_iterator(reinterpret_cast<CenterToOutIterator *>(
            memory + coefficient_memory_size_in_bytes)),
        partdist(reinterpret_cast<FL *>(
            memory + coefficient_memory_size_in_bytes + children_iterator_size_in_bytes)),
        center_partsums(reinterpret_cast<FL *>(
            memory + coefficient_memory_size_in_bytes + children_iterator_size_in_bytes + partdist_memory_size_in_bytes)),
        open_node_count(reinterpret_cast<unsigned int *>(
            memory + coefficient_memory_size_in_bytes + children_iterator_size_in_bytes +
            partdist_memory_size_in_bytes + center_partsums_memory_size_in_bytes))
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
     coefficient_memory_size_in_bytes + children_iterator_size_in_bytes + partdist_memory_size_in_bytes + center_partsums_memory_size_in_bytes + open_node_count_memory_size_in_bytes;

  static_assert(memory_size_in_bytes_per_group < std::numeric_limits<unsigned int>::max(),
                "Requires more memory than indexable with unsigned int");

  __device__ __host__ inline unsigned int get_node_count(unsigned int tree_level)
  {
    assert(tree_level < levels);
    return open_node_count[tree_level];
  }

  /**
   * Returns the squared norm of the node_index-th point at level tree_level projected into
   * the subspace perpendicular to the first levels - tree_level basis vectors
   */
  __device__ __host__ inline FL get_partdist(unsigned int tree_level,
                                                                unsigned int node_index)
  {
    assert(tree_level < levels);
    assert(node_index < max_nodes_per_level);
    return partdist[tree_level * max_nodes_per_level + node_index];
  }

  /**
   * Sets the squared norm of the node_index-th point at level tree_level projected into
   * the subspace perpendicular to the first levels - tree_level basis vectors
   */
  __device__ __host__ inline void set_partdist(unsigned int tree_level,
                                                                  unsigned int node_index, FL value)
  {
    assert(tree_level < levels);
    assert(node_index < max_nodes_per_level);
    partdist[tree_level * max_nodes_per_level + node_index] = value;
  }

  /**
   * Returns the inner product of the node_index-th point at level tree_level with
   * the given orthogonalized basis vector (given via index in basis)
   */
  __device__ __host__ inline FL get_center_partsum(unsigned int tree_level, unsigned int node_index,
                                                  unsigned int orthogonalized_basis_vector_index)
  {
    assert(tree_level < levels);
    assert(node_index < max_nodes_per_level);
    assert(orthogonalized_basis_vector_index < levels);
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
    assert(tree_level < levels);
    assert(node_index < max_nodes_per_level);
    assert(orthogonalized_basis_vector_index < levels);
    center_partsums[tree_level * levels * max_nodes_per_level +
                    orthogonalized_basis_vector_index * max_nodes_per_level + node_index] = value;
  }

  /**
   * Returns the iterator over the children coefficients
   */
  __device__ __host__ inline CenterToOutIterator &get_children_iterator(unsigned int tree_level,
                                                                        unsigned int node_index)
  {
    assert(tree_level < levels);
    assert(node_index < max_nodes_per_level);
    return children_iterator[tree_level * max_nodes_per_level + node_index];
  }

  /**
   * Returns the coefficient of the node_index-th point at level tree_level w.r.t.
   * the given basis vector (given by index). Called x in other fplll algorithms
   */
  __device__ __host__ inline FL get_coefficient(unsigned int tree_level, unsigned int node_index,
                                                unsigned int coordinate)
  {
    assert(tree_level < levels);
    assert(node_index < max_nodes_per_level);
    assert(coordinate < dimension);
    return coefficients[tree_level * max_nodes_per_level * levels +
                        coordinate * max_nodes_per_level + node_index];
  }
  
  /**
   * Sets the coefficient of the node_index-th point at topmost tree level (i.e. 0) w.r.t.
   * the given basis vector (given by index). Use only for coefficients determined by the
   * starting point (i.e. coordinates + 1 >= levels)
   */
  __device__ __host__ inline void set_start_coefficient(unsigned int node_index, unsigned int coordinate, FL value) 
  {
    assert(coordinate + 1 >= levels);
    assert(coordinate < dimension);
    assert(node_index < max_nodes_per_level);
    coefficients[0 + coordinate * max_nodes_per_level + node_index] = value;
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
    assert(tree_level < levels);

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
      for (unsigned int i = levels - tree_level; i < dimension; ++i)
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

  template <typename CG, unsigned int block_size>
  __device__ __host__ inline unsigned int
  add_starting_task(CG &cooperative_group, PrefixCounter<CG, block_size> &prefix_counter, bool this_thread_active)
  {
    unsigned int total_new_tasks      = 0;
    const unsigned int existing_tasks = get_node_count(0);
    const unsigned int new_task_offset =
        prefix_counter.prefix_count(cooperative_group, this_thread_active, total_new_tasks);
    assert(existing_tasks + total_new_tasks <= max_nodes_per_level);
    const unsigned int new_task_index = new_task_offset + existing_tasks;

    if (cooperative_group.thread_rank() == 0)
    {
      assert(total_new_tasks + open_node_count[0] <= max_nodes_per_level);
      open_node_count[0] += total_new_tasks;
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
    assert(tree_level < levels);
    assert(active_thread_count <= open_node_count[tree_level]);
    assert(node_index ==
           open_node_count[tree_level] - active_thread_count + cooperative_group.thread_rank());

    unsigned int kept_tasks = 0;
    const bool is_active =
        keep_this_thread_task && cooperative_group.thread_rank() < active_thread_count;
    const unsigned int new_offset =
        prefix_counter.prefix_count(cooperative_group, is_active, kept_tasks);
    const unsigned int new_index = new_offset + open_node_count[tree_level] - active_thread_count;

    const unsigned int old_node_index = tree_level * max_nodes_per_level + node_index;
    const unsigned int new_node_index = tree_level * max_nodes_per_level + new_index;

    FL projected_point_norm_square_tmp;
    FL partial_sum_tmp[levels];
    int coefficient_tmp[dimension];
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
      partdist[new_node_index]          = projected_point_norm_square_tmp;
      children_iterator[new_node_index] = children_iterator_tmp;
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