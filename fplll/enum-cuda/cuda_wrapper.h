#ifndef FPLLL_CUDA_ENUM
#define FPLLL_CUDA_ENUM

#include <utility>
#include <type_traits>
#include <functional>
#include "memory.h"

namespace cuenum
{

typedef std::function<float(double, double*)> process_sol_fn;

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

constexpr CudaEnumOpts default_opts = {50, .5, 3, 8, 32 * 256};

std::vector<uint64_t> search_enumeration_cuda(const double *mu, const double *rdiag,
                                 const unsigned int enum_dimensions,
                                 const double *start_point_coefficients, unsigned int start_point_count,
                                 unsigned int start_point_dim, process_sol_fn evaluator,
                                 double initial_radius, CudaEnumOpts opts = default_opts);

/**
 * Allocates memory and fills it with the given start points, so that it is directly copyable to the device
 * memory space (and therefore a correct parameter for search_enumeration_cuda()). The start points are given
 * as an iterator over indexable objects, each containing the start_point_dim coefficients of one start point.
 * The memory is allocated in page-locked memory to improve copy efficiency, but the provided unique pointer 
 * will correctly free it.
 */
template <typename InputIt>
inline PinnedPtr<double>
create_start_point_array(size_t start_point_count, size_t start_point_dim,
                         InputIt begin, InputIt end)
{
  PinnedPtr<double> result = allocatePinnedMemory<double>(start_point_count * start_point_dim);
  size_t i = 0;
  for (InputIt it = begin; it != end; ++it)
  {
    const auto& point = it->second;
    for (size_t j = 0; j < start_point_dim; ++j)
    {
      result.get()[i * start_point_dim + j] = static_cast<double>(point[j].get_d());
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

constexpr size_t cudaenum_return_array_size = 256;

typedef void(extenum_cb_set_config)(double *mu, size_t mudim, bool mutranspose, double *rdiag,
                                    double *pruning);
typedef double(extenum_cb_process_sol)(double dist, double *sol);
typedef void(extenum_cb_process_subsol)(double dist, double *subsol, int offset);

std::array<uint64_t, cudaenum_return_array_size> ext_cuda_enumerate(const int dim, double maxdist, std::function<extenum_cb_set_config> cbfunc,
    std::function<extenum_cb_process_sol> cbsol, std::function<extenum_cb_process_subsol> cbsubsol,
    bool dual = false, bool findsubsols = false);

#endif