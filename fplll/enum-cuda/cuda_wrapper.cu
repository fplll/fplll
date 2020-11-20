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
inline std::vector<uint64_t> search_enumeration_choose_levels(
    const enumf *mu, const enumf *rdiag, const unsigned int enum_levels,
    const enumi *start_point_coefficients, unsigned int start_point_count,
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
  return enumerate(mu, rdiag, start_point_coefficients, start_point_dim, start_point_count, initial_radius,
            evaluator, opts);
}

template <int dimensions_per_level, int min_levels, int delta_levels>
inline std::vector<uint64_t> search_enumeration_choose_levels(const enumf *mu, const enumf *rdiag,
                                             const unsigned int enum_levels,
                                             const enumi *start_point_coefficients,
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
    return search_enumeration_choose_levels(mu, rdiag, enum_levels, start_point_coefficients,
                                     start_point_count, start_point_dim, initial_radius, evaluator,
                                     enum_opts, int_marker<dimensions_per_level>(),
                                     int_marker<min_levels>(), int_marker<delta_mid>());
  }
  else
  {
    return search_enumeration_choose_levels(
        mu, rdiag, enum_levels, start_point_coefficients, start_point_count, start_point_dim,
        initial_radius, evaluator, enum_opts, int_marker<dimensions_per_level>(),
        int_marker<min_levels + delta_mid + 1>(), int_marker<delta_levels - delta_mid - 1>());
  }
}

template <int dimensions_per_level>
inline std::vector<uint64_t> search_enumeration_choose_dims_per_level(
    const enumf *mu, const enumf *rdiag, const unsigned int enum_dimensions,
    const enumi *start_point_coefficients, unsigned int start_point_count,
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

  return search_enumeration_choose_levels(mu, rdiag, enum_levels, start_point_coefficients,
                                   start_point_count, start_point_dim, initial_radius, evaluator,
                                   enum_opts, int_marker<dimensions_per_level>(), int_marker<1>(),
                                   int_marker<cudaenum_max_levels - 1>());
}

template <int min_dimensions_per_level, int delta_dimensions_per_level>
inline std::vector<uint64_t> search_enumeration_choose_dims_per_level(
    const enumf *mu, const enumf *rdiag, const unsigned int enum_dimensions,
    const enumi *start_point_coefficients, unsigned int start_point_count,
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
    return search_enumeration_choose_dims_per_level(
        mu, rdiag, enum_dimensions, start_point_coefficients, start_point_count, start_point_dim,
        initial_radius, evaluator, enum_opts, int_marker<min_dimensions_per_level>(),
        int_marker<delta_mid>());
  }
  else
  {
    return search_enumeration_choose_dims_per_level(
        mu, rdiag, enum_dimensions, start_point_coefficients, start_point_count, start_point_dim,
        initial_radius, evaluator, enum_opts, int_marker<min_dimensions_per_level + delta_mid + 1>(),
        int_marker<delta_dimensions_per_level - delta_mid - 1>());
  }
}

std::vector<uint64_t> search_enumeration_cuda(const double *mu, const double *rdiag,
                             const unsigned int enum_dimensions,
                             const enumi *start_point_coefficients, unsigned int start_point_count,
                             unsigned int start_point_dim, process_sol_fn evaluator,
                             double initial_radius, CudaEnumOpts opts)
{
  return search_enumeration_choose_dims_per_level(mu, rdiag, enum_dimensions, start_point_coefficients,
                                           start_point_count, start_point_dim, initial_radius,
                                           evaluator, opts, int_marker<1>(),
                                           int_marker<cudaenum_max_dims_per_level - 1>());
}

}

// TODO: replace with recursive enumeration from fplll
#include "util.h"

PinnedPtr<enumi> enumerate_start_points(const int dim, const int start_dims, double radius_squared, const enumf* mu, const enumf* rdiag, unsigned int& start_point_count) {

  // later, when using an evaluator + recursive fplll enumeration, this pair interface is what we need
  typedef int Dummy;
  std::vector<std::pair<Dummy, std::vector<FloatWrapper>>> start_points;

  std::unique_ptr<enumf[]> lattice(new enumf[dim * dim]);
  std::memcpy(lattice.get(), mu, dim * dim * sizeof(enumf));
  for (unsigned int i = 0; i < dim; ++i) {
    for (unsigned int j = 0; j < dim; ++j) {
      lattice[i * dim + j] *= std::sqrt(rdiag[i]);
    }
  }

  std::vector<FloatWrapper> x;
  x.resize(start_dims, { 0 });
  x[start_dims - 1] = -1;
  std::function<void(const std::vector<FloatWrapper> &)> callback =
      [&start_points](const std::vector<FloatWrapper> &a) { 
        start_points.push_back(std::make_pair(0, a)); 
      };
  do
  {
    ++x[start_dims - 1];
  } while (!naive_enum_recursive<FloatWrapper, enumf>(x, start_dims, dim, 0, 0, lattice.get(), dim - 1,
                                                       radius_squared, callback));
  start_point_count = start_points.size();
  return cuenum::create_start_point_array(start_points.size(), start_dims, start_points.begin(), start_points.end());
}

std::array<uint64_t, cudaenum_return_array_size> ext_cuda_enumerate(const int dim, enumf maxdist, std::function<extenum_cb_set_config> cbfunc,
  std::function<extenum_cb_process_sol> cbsol, std::function<extenum_cb_process_subsol> cbsubsol,
  bool dual, bool findsubsols) 
{
  if (dual) {
    throw "Unsupported operation: dual == true";
  } else if (findsubsols) {
    throw "Unsupported operation: findsubsols == true";
  }
  std::array<uint64_t, cudaenum_return_array_size> result = {};
  PinnedPtr<enumf> mu = allocatePinnedMemory<enumf>(dim * dim);
  PinnedPtr<enumf> rdiag = allocatePinnedMemory<enumf>(dim);
  enumf radius = std::sqrt(maxdist);

  cbfunc(mu.get(), dim, true, rdiag.get(), nullptr);

  cuenum::CudaEnumOpts opts = cuenum::default_opts;
  
  int start_dims = 6;
  while ((dim - start_dims) % opts.dimensions_per_level != 0) {
    ++start_dims;
  }
  if (start_dims >= dim) {
      // use fallback, as cuda enumeration in such small dimensions is too much overhead
      return result;
  }

  unsigned int start_point_count = 0;
  PinnedPtr<enumi> start_point_array = enumerate_start_points(dim, start_dims, maxdist, mu.get(), rdiag.get(), start_point_count);

  std::vector<uint64_t> node_counts = cuenum::search_enumeration_cuda(mu.get(), rdiag.get(), dim - start_dims, start_point_array.get(), 
    start_point_count, start_dims, cbsol, radius, opts);

  std::copy(node_counts.begin(), node_counts.end(), result.begin());
  return result;
}