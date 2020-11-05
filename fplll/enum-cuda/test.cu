#include "cuda_runtime.h"
#include <array>
#include <functional>
#include <vector>

#include "wrapper.h"
#include "atomic.h"
#include "cuda_util.cuh"
#include "enum.cuh"

#include "testdata.h"
#include "util.h"

typedef float enumi;
typedef double enumf;

struct Evaluator
{
  __device__ __host__ void operator()(int x, unsigned int i, bool done, double norm_squared,
                                      uint32_t *enum_bound_location)
  {
    // printf("%d, ", x);
    if (done)
    {
      // printf("\n");
      uint32_t squared_norm_repr = float_to_int_order_preserving_bijection(norm_squared);
      uint32_t old_repr          = atomic_min(enum_bound_location, squared_norm_repr);
    }
  }
};

template <typename eval_sol_fn, unsigned int total_dims, unsigned int start_point_dim>
void search_arr(const std::array<std::array<float, total_dims>, total_dims> &mu,
                const std::vector<std::array<enumi, start_point_dim>> &start_points,
                eval_sol_fn callback)
{
  constexpr unsigned int mu_n = total_dims;
  constexpr unsigned int dimensions_per_level = 4;
  static_assert((total_dims - start_point_dim) % dimensions_per_level == 0,
                "enumerated dims must be dividable by dimensions_per_level");
  assert(opts.dimensions_per_level == dimensions_per_level);

  std::unique_ptr<enumf[]> host_mu(new enumf[mu_n * mu_n]);
  std::unique_ptr<enumf[]> host_rdiag(new enumf[mu_n]);
  for (unsigned int i = 0; i < mu_n; ++i)
  {
    host_rdiag[i] = mu[i][i];
    for (unsigned int j = 0; j < mu_n; ++j)
    {
      host_mu[i * mu_n + j] = mu[i][j] / host_rdiag[i];
    }
    host_rdiag[i] = host_rdiag[i] * host_rdiag[i];
  }
  CudaEnumOpts opts = default_opts;
  opts.dimensions_per_level = dimensions_per_level;
  cudaenum::Opts<(total_dims - start_point_dim) / dimensions_per_level, dimensions_per_level, 2048>
      static_opts = {
      opts.max_subtree_paths, opts.min_active_parents_percentage,
      2048 - 32 * opts.max_subtree_paths, opts.initial_nodes_per_group, opts.thread_count};

  cudaenum::enumerate(host_mu.get(), host_rdiag.get(), start_points[0].data(), start_point_dim,
                      start_points.size(), find_initial_radius(mu), callback, static_opts);
}

inline void gpu_test()
{
  constexpr unsigned int total_dim       = 50;
  constexpr unsigned int start_point_dim = 6;

  const std::array<std::array<float, total_dim>, total_dim> &mu = test_mu_knapsack_big;
  std::vector<std::array<enumi, start_point_dim>> start_points;

  std::array<enumi, start_point_dim> x;
  x[start_point_dim - 1] = -1;
  enumf radius           = find_initial_radius(mu) * 1.5;
  std::function<void(const std::array<float, start_point_dim> &)> callback =
      [&start_points](const auto &a) { start_points.push_back(a); };
  do
  {
    ++x[start_point_dim - 1];
  } while (!naive_enum_recursive<start_point_dim, total_dim>(x, 0, 0, mu, total_dim - 1,
                                                             radius * radius, callback));

  Evaluator evaluator;
  search_arr(mu, start_points, evaluator);
}

// inline void cpu_test4d()
//{
//  constexpr unsigned int levels               = 1;
//  constexpr unsigned int dimensions_per_level = 4;
//  constexpr unsigned int dimensions           = levels * dimensions_per_level;
//  constexpr unsigned int max_nodes_per_level  = 100;
//
//  std::array<std::array<float, dimensions>, dimensions> host_mu = {
//      {{3, 2, 4, 3}, {0, 3, 4, 2}, {0, 0, 3, 3}, {0, 0, 0, 2}}};
//
//  single_thread group;
//  PrefixCounter<single_thread, 1> prefix_counter;
//
//  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> buffer(
//      new unsigned char[SubtreeEnumerationBuffer<levels, dimensions_per_level,
//                                                 max_nodes_per_level>::memory_size_in_bytes]);
//  buffer.init(group);
//
//  std::unique_ptr<enumf[]> mu(new enumf[dimensions * dimensions]);
//  std::unique_ptr<enumf[]> local_mu(new enumf[dimensions_per_level * dimensions_per_level]);
//  std::unique_ptr<enumf[]> rdiag(new enumf[dimensions]);
//  for (unsigned int i = 0; i < dimensions; ++i)
//  {
//    rdiag[i] = host_mu[i][i];
//    for (unsigned int j = 0; j < dimensions; ++j)
//    {
//      mu[i * dimensions + j] = host_mu[i][j] / rdiag[i];
//    }
//    rdiag[i] = rdiag[i] * rdiag[i];
//  }
//  enumf radius             = 3.01;
//  uint32_t radius_location = float_to_int_order_preserving_bijection(radius * radius);
//
//  const unsigned int index = buffer.add_node(0, 0);
//  for (unsigned int i = 0; i < dimensions; ++i)
//  {
//    buffer.set_center_partsum(0, index, i, 0);
//  }
//  buffer.init_subtree(0, index, 0, 0);
//
//  unsigned int counter            = 0;
//  unsigned long long node_counter = 0;
//  StrategyOpts opts               = {5, .5, max_nodes_per_level / 2};
//  clear_level(
//      group, prefix_counter, &counter, buffer, 0, Matrix(mu.get(), dimensions), rdiag.get(),
//      &radius_location, [](...) {}, opts, PerfCounter(&node_counter));
//}

int main()
{
  gpu_test();

  return 0;
}
