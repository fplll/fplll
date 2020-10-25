#include "enum.cu"
#include "enum2.cu"
#include "enum.cuh"
#include "testdata.h"

#include <array>
#include <functional>

void simple_gpu_test() {

  constexpr unsigned int dimensions         = 3;
  constexpr unsigned int levels             = dimensions;

  std::array<std::array<float, dimensions>, dimensions> host_mu = {
      {{4.58257569, 4.3643578, 3.27326835}, {0., 1.71824939, 2.16166858}, {0., 0., 1.27000127}}};

  std::vector<std::array<float, dimensions - levels + 1>> start_points;
  start_points.push_back({0});
  start_points.push_back({1});

  unsigned int output_point_count = 0;
  CudaPtr<float> shortest_points_per_thread = search_enumeration<float, levels, dimensions>(
      host_mu, start_points, output_point_count, find_initial_radius(host_mu) * 1.1);

  std::unique_ptr<float[]> result(new float[dimensions * output_point_count]);
  check(cudaMemcpy(result.get(), shortest_points_per_thread.get(), levels * output_point_count * sizeof(float), cudaMemcpyDeviceToHost));

  std::array<float, levels> expected = { 1, -1, 0 };
  if (!contains_solution(result.get(), expected.data(), levels, output_point_count)) {
      throw "Test failed";
  }
}

void complex_gpu_test() {

  constexpr unsigned int dimensions = 50;
  constexpr unsigned int levels     = dimensions;

  std::vector<std::array<float, dimensions - levels + 1>> start_points;
  start_points.push_back({0});
  start_points.push_back({1});

  unsigned int output_point_count = 0;
  CudaPtr<float> shortest_points_per_thread = search_enumeration<float, levels, dimensions>(
      test_mu_normal, start_points, output_point_count,
      find_initial_radius(test_mu_normal) * 1.1);

  std::unique_ptr<float[]> result(new float[dimensions * output_point_count]);
  check(cudaMemcpy(result.get(), shortest_points_per_thread.get(), levels * output_point_count * sizeof(float), cudaMemcpyDeviceToHost));

  if (!contains_solution(result.get(), test_solution_normal.data(), levels, output_point_count)) {
      throw "Test failed";
  }
}

void complex_multiblock_gpu_test()
{
  constexpr unsigned int dimensions = 50;
  constexpr unsigned int levels     = dimensions - 2;

  std::vector<std::array<float, dimensions - levels + 1>> start_points;

  float initial_radius = find_initial_radius(test_mu_normal) * 1.1;

  std::array<float, 3> x;
  std::function<void(const std::array<float, 3> &)> fn = [&start_points](auto point) {
    start_points.emplace_back(point);
  };
  CenterToOutIterator iter;
  do
  {
    x[2] = *iter;
    ++iter;
  } while (!naive_enum_recursive<3, dimensions>(x, 0.f, 0.f, test_mu_normal, dimensions - 1,
                                                initial_radius * initial_radius, fn));

  unsigned int output_point_count           = 0;
  CudaPtr<float> shortest_points_per_thread = search_enumeration<float, levels, dimensions>(
      test_mu_normal, start_points, output_point_count, initial_radius);

  std::unique_ptr<float[]> result(new float[dimensions * output_point_count]);
  check(cudaMemcpy(result.get(), shortest_points_per_thread.get(),
                   levels * output_point_count * sizeof(float), cudaMemcpyDeviceToHost));

  if (!contains_solution(result.get(), test_solution_normal.data(), levels, output_point_count))
  {
    throw "Test failed";
  }
}

int main()
{
  hybrid_enum::test();
  /*simple_gpu_test();
  complex_multiblock_gpu_test();
  complex_gpu_test();*/

  return 0;
}
