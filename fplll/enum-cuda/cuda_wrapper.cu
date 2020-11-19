#include "atomic.h"
#include "enum.cuh"
#include "cuda_wrapper.h"

namespace cuenum {

constexpr int cudaenum_max_dims_per_level  = 4;
constexpr int cudaenum_max_levels          = 19;
constexpr unsigned int cudaenum_max_nodes_per_level = 3100;

template <int min> struct int_marker
{
};

template <int dimensions_per_level, int levels>
inline void search_enumeration_choose_levels(
    const enumf *mu, const enumf *rdiag, const unsigned int enum_levels,
    const float *start_point_coefficients, unsigned int start_point_count,
    unsigned int start_point_dim, enumf initial_radius, process_sol_fn evaluator,
    CudaEnumOpts enum_opts, int_marker<dimensions_per_level>, int_marker<levels>, int_marker<0>)
{
  if (enum_levels != levels) {
    throw "enumeration dimension must be within the allowed interval";
  }
  assert(enum_opts.max_subtree_paths * enumerate_cooperative_group_size <
         cudaenum_max_nodes_per_level);
  unsigned int max_children_node_count =
      cudaenum_max_nodes_per_level - enum_opts.max_subtree_paths * enumerate_cooperative_group_size;
  Opts<levels, dimensions_per_level, cudaenum_max_nodes_per_level> opts = {
      enum_opts.max_subtree_paths, enum_opts.min_active_parents_percentage, max_children_node_count,
      enum_opts.initial_nodes_per_group, enum_opts.thread_count};
  enumerate(mu, rdiag, start_point_coefficients, start_point_dim, start_point_count, initial_radius,
            evaluator, opts);
}

template <int dimensions_per_level, int min_levels, int delta_levels>
inline void search_enumeration_choose_levels(const enumf *mu, const enumf *rdiag,
                                             const unsigned int enum_levels,
                                             const float *start_point_coefficients,
                                             unsigned int start_point_count,
                                             unsigned int start_point_dim, enumf initial_radius,
                                             process_sol_fn evaluator, CudaEnumOpts enum_opts,
                                             int_marker<dimensions_per_level>,
                                             int_marker<min_levels>, int_marker<delta_levels>)
{
  static_assert(delta_levels >= 0, "delta_levels >= 0 must hold");
  assert(enum_levels >= min_levels);
  assert(enum_levels <= min_levels + delta_levels);

  constexpr unsigned int delta_mid = delta_levels / 2;
  if (enum_levels <= min_levels + delta_mid)
  {
    search_enumeration_choose_levels(mu, rdiag, enum_levels, start_point_coefficients,
                                     start_point_count, start_point_dim, initial_radius, evaluator,
                                     enum_opts, int_marker<dimensions_per_level>(),
                                     int_marker<min_levels>(), int_marker<delta_mid>());
  }
  else
  {
    search_enumeration_choose_levels(
        mu, rdiag, enum_levels, start_point_coefficients, start_point_count, start_point_dim,
        initial_radius, evaluator, enum_opts, int_marker<dimensions_per_level>(),
        int_marker<min_levels + delta_mid + 1>(), int_marker<delta_levels - delta_mid - 1>());
  }
}

template <int dimensions_per_level>
inline void search_enumeration_choose_dims_per_level(
    const enumf *mu, const enumf *rdiag, const unsigned int enum_dimensions,
    const float *start_point_coefficients, unsigned int start_point_count,
    unsigned int start_point_dim, enumf initial_radius, process_sol_fn evaluator,
    CudaEnumOpts enum_opts, int_marker<dimensions_per_level>, int_marker<0>)
{
  if (enum_opts.dimensions_per_level != dimensions_per_level) {
    throw "enum_opts.dimensions_per_level not within allowed interval";
  }
  if (enum_dimensions % dimensions_per_level != 0) {
    throw "enumeration dimension count (i.e. dimensions minus start point dimensions) must be divisible by enum_opts.dimensions_per_level";
  }
  unsigned int enum_levels = enum_dimensions / dimensions_per_level;

  search_enumeration_choose_levels(mu, rdiag, enum_levels, start_point_coefficients,
                                   start_point_count, start_point_dim, initial_radius, evaluator,
                                   enum_opts, int_marker<dimensions_per_level>(), int_marker<1>(),
                                   int_marker<cudaenum_max_levels - 1>());
}

template <int min_dimensions_per_level, int delta_dimensions_per_level>
inline void search_enumeration_choose_dims_per_level(
    const enumf *mu, const enumf *rdiag, const unsigned int enum_dimensions,
    const float *start_point_coefficients, unsigned int start_point_count,
    unsigned int start_point_dim, enumf initial_radius, process_sol_fn evaluator,
    CudaEnumOpts enum_opts, int_marker<min_dimensions_per_level>,
    int_marker<delta_dimensions_per_level>)
{
  static_assert(delta_dimensions_per_level >= 0, "delta_dimensions_per_level >= 0must hold");
  assert(enum_opts.dimensions_per_level >= min_dimensions_per_level);
  assert(enum_opts.dimensions_per_level <= min_dimensions_per_level + delta_dimensions_per_level);

  constexpr unsigned int delta_mid = delta_dimensions_per_level / 2;
  if (enum_opts.dimensions_per_level <= min_dimensions_per_level + delta_mid)
  {
    search_enumeration_choose_dims_per_level(
        mu, rdiag, enum_dimensions, start_point_coefficients, start_point_count, start_point_dim,
        initial_radius, evaluator, enum_opts, int_marker<min_dimensions_per_level>(),
        int_marker<delta_mid>());
  }
  else
  {
    search_enumeration_choose_dims_per_level(
        mu, rdiag, enum_dimensions, start_point_coefficients, start_point_count, start_point_dim,
        initial_radius, evaluator, enum_opts, int_marker<min_dimensions_per_level + delta_mid + 1>(),
        int_marker<delta_dimensions_per_level - delta_mid - 1>());
  }
}

void search_enumeration_cuda(const double *mu, const double *rdiag,
                             const unsigned int enum_dimensions,
                             const float *start_point_coefficients, unsigned int start_point_count,
                             unsigned int start_point_dim, process_sol_fn evaluator,
                             double initial_radius, CudaEnumOpts opts)
{
  search_enumeration_choose_dims_per_level(mu, rdiag, enum_dimensions, start_point_coefficients,
                                           start_point_count, start_point_dim, initial_radius,
                                           evaluator, opts, int_marker<1>(),
                                           int_marker<cudaenum_max_dims_per_level - 1>());
}

}