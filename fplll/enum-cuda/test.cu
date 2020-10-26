#include "enum.cu"
#include "testdata.h"

#include <array>
#include <functional>

inline void gpu_test()
{
  constexpr unsigned int levels               = 13;
  constexpr unsigned int dimensions_per_level = 4;
  constexpr unsigned int dimensions           = levels * dimensions_per_level;
  constexpr unsigned int max_nodes_per_level  = 800;
  constexpr unsigned int start_point_dim      = 8;
  constexpr unsigned int mu_n                 = dimensions + start_point_dim;

  const std::array<std::array<float, mu_n>, mu_n> &mu = test_mu_big;

  std::vector<std::array<enumi, start_point_dim>> start_points;
  std::array<enumi, start_point_dim> x;
  x[start_point_dim - 1] = -1;
  enumf radius           = find_initial_radius(mu) * 1.5;
  std::function<void(const std::array<float, start_point_dim> &)> callback =
      [&start_points](const auto &a) { start_points.push_back(a); };
  do
  {
    ++x[start_point_dim - 1];
  } while (!naive_enum_recursive<start_point_dim, mu_n>(x, 0, 0, mu, mu_n - 1,
                                                        radius * radius, callback));

  search<levels, dimensions_per_level, max_nodes_per_level>(mu, start_points, 20, 32, 8);
}

template<typename Dummy>
inline void gpu_test4d()
{
  constexpr unsigned int levels               = 3;
  constexpr unsigned int dimensions_per_level = 1;
  constexpr unsigned int dimensions           = levels * dimensions_per_level;
  constexpr unsigned int max_nodes_per_level  = 5000;
  constexpr unsigned int start_point_dim      = 1;
  constexpr unsigned int mu_n                 = dimensions + start_point_dim;

  std::array<std::array<float, mu_n>, mu_n> test_mu_tiny = {
      {{3, 2, 4, 3}, {0, 3, 4, 2}, {0, 0, 3, 3}, {0, 0, 0, 2}}};
  const std::array<std::array<float, mu_n>, mu_n> &mu = test_mu_tiny;

  std::vector<std::array<enumi, start_point_dim>> start_points;
  start_points.push_back({0});
  start_points.push_back({1});

  search<levels, dimensions_per_level, max_nodes_per_level>(mu, start_points, 10, 1, 1);
}

inline void cpu_test()
{
  constexpr unsigned int levels               = 5;
  constexpr unsigned int dimensions_per_level = 4;
  constexpr unsigned int dimensions           = levels * dimensions_per_level;
  constexpr unsigned int max_nodes_per_level  = 200000;

  const std::array<std::array<float, dimensions>, dimensions> &host_mu = test_mu_small;

  single_thread group;
  PrefixCounter<single_thread, 1> prefix_counter;

  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> buffer(
      new unsigned char[SubtreeEnumerationBuffer<levels, dimensions_per_level,
                                                 max_nodes_per_level>::memory_size_in_bytes]);
  buffer.init(group);

  std::unique_ptr<enumf[]> mu(new enumf[dimensions * dimensions]);
  std::unique_ptr<enumf[]> local_mu(new enumf[dimensions_per_level * dimensions_per_level]);
  std::unique_ptr<enumf[]> rdiag(new enumf[dimensions]);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    rdiag[i] = host_mu[i][i];
    for (unsigned int j = 0; j < dimensions; ++j)
    {
      mu[i * dimensions + j] = host_mu[i][j] / rdiag[i];
    }
    rdiag[i] = rdiag[i] * rdiag[i];
  }
  enumf radius             = find_initial_radius(host_mu) * 1.1;
  uint32_t radius_location = float_to_int_order_preserving_bijection(radius * radius);

  const unsigned int index = buffer.add_subtree(0, 0);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    buffer.set_center_partsum(0, index, i, 0);
  }
  buffer.init_subtree(0, index, 0, 0);

  unsigned int counter = 0;
  unsigned long long node_counter = 0;
  clear_level<single_thread, 1, levels, dimensions_per_level, max_nodes_per_level>(
      group, prefix_counter, reinterpret_cast<unsigned char *>(local_mu.get()), &counter, buffer, 0,
      mu.get(), dimensions, rdiag.get(), &radius_location, 5,
      PerfCounter(&node_counter));
}

inline void cpu_test4d()
{
  constexpr unsigned int levels               = 1;
  constexpr unsigned int dimensions_per_level = 4;
  constexpr unsigned int dimensions           = levels * dimensions_per_level;
  constexpr unsigned int max_nodes_per_level  = 100;

  std::array<std::array<float, dimensions>, dimensions> host_mu = {
      {{3, 2, 4, 3}, {0, 3, 4, 2}, {0, 0, 3, 3}, {0, 0, 0, 2}}};

  single_thread group;
  PrefixCounter<single_thread, 1> prefix_counter;

  SubtreeEnumerationBuffer<levels, dimensions_per_level, max_nodes_per_level> buffer(
      new unsigned char[SubtreeEnumerationBuffer<levels, dimensions_per_level,
                                                 max_nodes_per_level>::memory_size_in_bytes]);
  buffer.init(group);

  std::unique_ptr<enumf[]> mu(new enumf[dimensions * dimensions]);
  std::unique_ptr<enumf[]> local_mu(new enumf[dimensions_per_level * dimensions_per_level]);
  std::unique_ptr<enumf[]> rdiag(new enumf[dimensions]);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    rdiag[i] = host_mu[i][i];
    for (unsigned int j = 0; j < dimensions; ++j)
    {
      mu[i * dimensions + j] = host_mu[i][j] / rdiag[i];
    }
    rdiag[i] = rdiag[i] * rdiag[i];
  }
  enumf radius             = 3.01;
  uint32_t radius_location = float_to_int_order_preserving_bijection(radius * radius);

  const unsigned int index = buffer.add_subtree(0, 0);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    buffer.set_center_partsum(0, index, i, 0);
  }
  buffer.init_subtree(0, index, 0, 0);

  unsigned int counter      = 0;
  unsigned long long node_counter = 0;
  clear_level<single_thread, 1, levels, dimensions_per_level, max_nodes_per_level>(
      group, prefix_counter, reinterpret_cast<unsigned char *>(local_mu.get()), &counter, buffer, 0,
      mu.get(), dimensions, rdiag.get(), &radius_location, 5,
      PerfCounter(&node_counter));
}

int main()
{
  gpu_test();

  return 0;
}
