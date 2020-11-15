#ifndef FPLLL_CUDA_ENUM
#define FPLLL_CUDA_ENUM

#include <utility>
#include <type_traits>
#include "cuda_check.h"

namespace cuenum
{

enum struct EvaluatorStrategy
{
  FAST
};

struct CudaEnumOpts
{
  // maximal amount of paths that will be searched using recursive enumeration during each algorithm
  // step by each thread. When this is exceeded, the recursive enumeration state is stored and
  // resumed in the next step. Use for load balancing to prevent threads with small subtrees to wait
  // for threads with very big subtrees.
  unsigned int max_subtree_paths;
  // stop children generation when the percentage of parent points that have still unprocessed
  // children drops beneath this percentage
  float min_active_parents_percentage;
  // height of the subtrees that are searched using recursive enumeration.
  unsigned int dimensions_per_level;
  // how many start points should be assigned to each cooperative group
  unsigned int initial_nodes_per_group;
  // how many cuda threads to use for the search. If this is not a multiple of the block
  // size, it will be rounded up to one
  unsigned int thread_count;
};

constexpr CudaEnumOpts default_opts = {50, .5, 4, 8, 32 * 256};

void search_enumeration_cuda(const double *mu, const double *rdiag,
                             const unsigned int enum_dimensions,
                             const float *start_point_coefficients, unsigned int start_point_count,
                             unsigned int start_point_dim, EvaluatorStrategy evaluator,
                             double initial_radius, CudaEnumOpts opts = default_opts);

template <typename T> struct PinnedDeleter
{
  inline void operator()(T *ptr) const { check(cudaFreeHost(ptr)); }
};

/**
 * Allocates memory and fills it with the given start points, so that it is directly copyable to the device
 * memory space (and therefore a correct parameter for search_enumeration_cuda()). The start points are given
 * as an iterator over indexable objects, each containing the start_point_dim coefficients of one start point.
 * The memory is allocated in page-locked memory to improve copy efficiency, but the provided unique pointer 
 * will correctly free it.
 */
template <typename InputIt>
inline std::unique_ptr<float, PinnedDeleter<float>>
create_start_point_array(size_t start_point_count, size_t start_point_dim,
                         InputIt begin, InputIt end)
{
  float *ptr;
  check(cudaMallocHost(&ptr, start_point_count * start_point_dim * sizeof(float)));
  std::unique_ptr<float, PinnedDeleter<float>> result(ptr, PinnedDeleter<float>());
  size_t i = 0;
  for (InputIt it = begin; it != end; ++it)
  {
    const auto& point = it->second;
    for (size_t j = 0; j < start_point_dim; ++j)
    {
      result.get()[i * start_point_dim + j] = static_cast<float>(point[j].get_d());
    }
    ++i;
  }
  if (i != start_point_count)
  {
    throw "Given and actual start point counts do not match";
  }
  return result;
}

}  // namespace cuenum

#endif