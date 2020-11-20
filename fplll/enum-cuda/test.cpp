#include <array>
#include <functional>
#include <vector>
#include <memory>

#include "cuda_wrapper.h"
#include "atomic.h"

#include "testdata.h"
#include "util.h"

using namespace cuenum;

typedef float enumi;
typedef double enumf;

template <unsigned int total_dims, unsigned int start_point_dim>
void search_arr_dyn(const std::array<std::array<float, total_dims>, total_dims> &mu,
                    const std::vector<std::pair<int, std::vector<FloatWrapper>>> &start_points)
{
  constexpr unsigned int mu_n                 = total_dims;
  constexpr unsigned int dimensions_per_level = 3;
  static_assert((total_dims - start_point_dim) % dimensions_per_level == 0,
                "enumerated dims must be dividable by dimensions_per_level");

  std::unique_ptr<enumf[]> host_mu(new enumf[mu_n * mu_n]);
  std::unique_ptr<enumf[]> host_rdiag(new enumf[mu_n]);
  for (size_t i = 0; i < mu_n; ++i)
  {
    host_rdiag[i] = mu[i][i];
    for (size_t j = 0; j < mu_n; ++j)
    {
      host_mu[i * mu_n + j] = mu[i][j] / host_rdiag[i];
    }
    host_rdiag[i] = host_rdiag[i] * host_rdiag[i];
  }

  CudaEnumOpts opts         = default_opts;
  opts.dimensions_per_level = dimensions_per_level;

  process_sol_fn callback = [](double norm, float* x) -> float { return norm; };

  auto start_point_memory = create_start_point_array(start_points.size(), start_point_dim,
                                                     start_points.begin(), start_points.end());
  search_enumeration_cuda(host_mu.get(), host_rdiag.get(), total_dims - start_point_dim,
                          start_point_memory.get(), static_cast<unsigned int>(start_points.size()),
                          start_point_dim, callback,
                          find_initial_radius<float, total_dims>(mu) * 1.1, opts);
}

void test_direct()
{
  constexpr unsigned int total_dim       = 50;
  constexpr unsigned int start_point_dim = 8;

  const std::array<std::array<float, total_dim>, total_dim> &mu = test_mu_knapsack_big;
  std::vector<std::pair<int, std::vector<FloatWrapper>>> start_points;

  std::vector<FloatWrapper> x;
  x.resize(start_point_dim, { 0 });
  x[start_point_dim - 1] = -1;
  enumf radius           = find_initial_radius<float, total_dim>(mu) * 1.1;
  std::function<void(const std::vector<FloatWrapper> &)> callback =
      [&start_points](const std::vector<FloatWrapper> &a) { start_points.push_back(std::make_pair(0, a)); };
  do
  {
    ++x[start_point_dim - 1];
  } while (!naive_enum_recursive<FloatWrapper, float>(x, start_point_dim, total_dim, 0, 0, &mu[0][0], total_dim - 1,
                                                             static_cast<float>(radius * radius), callback));
  
  search_arr_dyn<total_dim, start_point_dim>(mu, start_points);
}

void test_fplll_like() {
  
  constexpr unsigned int total_dim       = 50;
  const std::array<std::array<float, total_dim>, total_dim> &lattice = test_mu_knapsack_big;

  double maxdist = find_initial_radius<float, total_dim>(lattice) * 1.1;
  maxdist = maxdist * maxdist; 
  std::function<extenum_cb_set_config> set_config = [&lattice](double *mu, size_t mudim, bool mutranspose, double *rdiag, double *pruning) {
    assert(mutranspose);
    for (unsigned int i = 0; i < mudim; ++i) {
      for (unsigned int j = i; j < mudim; ++j) {
        mu[i * mudim + j] = lattice[i][j] / lattice[i][i];
      }
      rdiag[i] = lattice[i][i] * lattice[i][i];
    }
  };
  std::function<extenum_cb_process_sol> process_sol = [](double norm_square, float* x)-> double { return norm_square; };
  ext_cuda_enumerate(total_dim, maxdist, set_config, process_sol, nullptr);
}

int main()
{
  test_fplll_like();

  return 0;
}
