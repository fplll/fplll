#include "wrapper.h"
#include "atomic.h"
#include "enum.cuh"

constexpr unsigned int cudaenum_max_dims_per_level   = 5;
constexpr unsigned int cudaenum_max_levels           = 15;
constexpr unsigned int cudaenum_max_nodes_per_level  = 3100;

template <typename eval_sol_fn, unsigned int min_level, unsigned int max_level, unsigned int dimensions_per_level>
inline void search_enumeration_choose_levels(const enumf *mu, const enumf *rdiag,
                                             const unsigned int enum_levels,
                                             const float *start_point_coefficients,
                                             unsigned int start_point_count,
                                             unsigned int start_point_dim, enumf initial_radius,
                                             eval_sol_fn callback, CudaEnumOpts enum_opts)
{
  static_assert(min_level <= max_level, "min_level <= max_level must hold");
  assert(enum_levels >= min_level);
  assert(enum_levels <= max_level);
  if constexpr (min_level == max_level)
  {
    assert(enum_opts.max_subtree_paths * cudaenum::enumerate_cooperative_group_size <
           cudaenum_max_nodes_per_level);
    unsigned int max_children_node_count =
        cudaenum_max_nodes_per_level -
        enum_opts.max_subtree_paths * cudaenum::enumerate_cooperative_group_size;
    cudaenum::Opts<min_level, dimensions_per_level, cudaenum_max_nodes_per_level> opts = {
        enum_opts.max_subtree_paths, enum_opts.min_active_parents_percentage,
        max_children_node_count, enum_opts.initial_nodes_per_group, enum_opts.thread_count};
    cudaenum::enumerate(mu, rdiag, start_point_coefficients, start_point_dim, start_point_count,
                        initial_radius, callback, opts);
  }
  else
  {
    constexpr unsigned int mid = (min_level + max_level) / 2;
    if (enum_levels <= mid)
    {
      search_enumeration_choose_levels<eval_sol_fn, min_level, mid, dimensions_per_level>(
          mu, rdiag, enum_levels, start_point_coefficients, start_point_count, start_point_dim,
          initial_radius, callback, enum_opts);
    }
    else
    {
      search_enumeration_choose_levels<eval_sol_fn, mid + 1, max_level, dimensions_per_level>(
          mu, rdiag, enum_levels, start_point_coefficients, start_point_count, start_point_dim,
          initial_radius, callback, enum_opts);
    }
  }
}


template <typename eval_sol_fn, unsigned int min_dimensions_per_level, unsigned int max_dimensions_per_level>
inline void search_enumeration_choose_dims_per_level(const enumf *mu, const enumf *rdiag,
                                                     const unsigned int enum_levels,
                                                     const float *start_point_coefficients,
                                                     unsigned int start_point_count,
                                                     unsigned int start_point_dim,
                                                     enumf initial_radius, eval_sol_fn callback,
                                                     CudaEnumOpts enum_opts)
{
  static_assert(min_dimensions_per_level <= max_dimensions_per_level,
                "min_dimensions_per_level <= max_dimensions_per_level must hold");
  assert(enum_opts.dimensions_per_level >= min_dimensions_per_level);
  assert(enum_opts.dimensions_per_level <= max_dimensions_per_level);
  if constexpr (min_dimensions_per_level == max_dimensions_per_level)
  {
    search_enumeration_choose_levels<eval_sol_fn, 1, cudaenum_max_levels, min_dimensions_per_level>(
        mu, rdiag, enum_levels, start_point_coefficients, start_point_count, start_point_dim,
        initial_radius, callback, enum_opts);
  }
  else
  {
    constexpr unsigned int mid = (min_dimensions_per_level + max_dimensions_per_level) / 2;
    if (enum_opts.dimensions_per_level <= mid)
    {
      search_enumeration_choose_dims_per_level<eval_sol_fn, min_dimensions_per_level, mid>(
          mu, rdiag, enum_levels, start_point_coefficients, start_point_count, start_point_dim,
          initial_radius, callback, enum_opts);
    }
    else
    {
      search_enumeration_choose_dims_per_level<eval_sol_fn, mid + 1, max_dimensions_per_level>(
          mu, rdiag, enum_levels, start_point_coefficients, start_point_count, start_point_dim,
          initial_radius, callback, enum_opts);
    }
  }
}

struct FastEvaluator
{
  enumf *norm_squares;
  enumi *result;
  unsigned int thread_count;

  inline FastEvaluator(enumf *norm_squared, enumi *result, unsigned int thread_count)
      : norm_squares(norm_squares), result(result), thread_count(thread_count)
  {
  }

  __device__ __host__ inline void operator()(enumi x, unsigned int i, bool done, enumf norm_square,
                                             uint32_t *enum_bound_location)
  {
    if (i == 0)
    {
      uint32_t squared_norm_repr = float_to_int_order_preserving_bijection(norm_square);
      uint32_t old_repr          = atomic_min(enum_bound_location, squared_norm_repr);
      norm_squares[thread_id()]  = norm_square;
    }
    result[i * thread_count + thread_id()] = x;
  }
};

inline void search_enumeration_choose_evaluator(const enumf *mu, const enumf *rdiag,
                                                const unsigned int enum_levels,
                                                const float *start_point_coefficients,
                                                unsigned int start_point_count,
                                                unsigned int start_point_dim,
                                                EvaluatorStrategy evaluator, enumf initial_radius,
                                                CudaEnumOpts opts)
{
  const unsigned int thread_count = cudaenum::get_started_thread_count(opts.thread_count);
  const unsigned int total_dimension = start_point_dim + enum_levels * opts.dimensions_per_level;
  switch (evaluator)
  {
  case EvaluatorStrategy::FAST:
    typedef FastEvaluator eval_sol_fn;
    CudaPtr<enumf> norms_squared = alloc(enumf, thread_count);
    CudaPtr<enumi> coefficients  = alloc(enumi, thread_count * total_dimension);
    FastEvaluator evaluator(norms_squared.get(), coefficients.get(), thread_count);
    search_enumeration_choose_dims_per_level<eval_sol_fn, 1, cudaenum_max_dims_per_level>(
        mu, rdiag, enum_levels, start_point_coefficients, start_point_count, start_point_dim,
        initial_radius, evaluator, opts);
    break;
  }
}

void search_enumeration_cuda(const enumf *mu, const enumf *rdiag,
                             const unsigned int enum_dimensions,
                             const float *start_point_coefficients, unsigned int start_point_count,
                             unsigned int start_point_dim, EvaluatorStrategy evaluator,
                             enumf initial_radius, CudaEnumOpts opts)
{
  assert(enum_dimensions % opts.dimensions_per_level == 0);
  unsigned int enum_levels = enum_dimensions / opts.dimensions_per_level;

  search_enumeration_choose_evaluator(mu, rdiag, enum_levels, start_point_coefficients,
                                      start_point_count, start_point_dim, evaluator, initial_radius,
                                      opts);

}