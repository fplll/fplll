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

struct FT {
  float value;

  inline double get_d() const {
    return value;
  }

  inline FT operator-(float rhs) const {
    return { value - rhs };
  }

  inline FT operator*(float rhs) const {
    return { value * rhs };
  }

  inline FT& operator=(float rhs) {
    value = rhs;
    return *this;
  }

  inline FT& operator++() {
    ++value;
    return *this;
  }

  inline operator float() const {
    return value;
  }
};

template <unsigned int total_dims, unsigned int start_point_dim>
void search_arr_dyn(const std::array<std::array<float, total_dims>, total_dims> &mu,
                    const std::vector<std::pair<int, std::array<FT, start_point_dim>>> &start_points)
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

inline void gpu_test()
{
  constexpr unsigned int total_dim       = 60;
  constexpr unsigned int start_point_dim = 6;

  const std::array<std::array<float, total_dim>, total_dim> &mu = test_mu_big;
  std::vector<std::pair<int, std::array<FT, start_point_dim>>> start_points;

  std::array<FT, start_point_dim> x;
  x[start_point_dim - 1] = -1;
  enumf radius           = find_initial_radius<float, total_dim>(mu) * 1.5;
  std::function<void(const std::array<FT, start_point_dim> &)> callback =
      [&start_points](const std::array<FT, start_point_dim> &a) { start_points.push_back(std::make_pair(0, a)); };
  do
  {
    ++x[start_point_dim - 1];
  } while (!naive_enum_recursive<FT, start_point_dim, total_dim>(x, 0, 0, mu, total_dim - 1,
                                                             static_cast<float>(radius * radius), callback));
  
  search_arr_dyn<total_dim, start_point_dim>(mu, start_points);
}

int main()
{
  gpu_test();

  return 0;
}
