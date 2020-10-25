#pragma once
#include "cooperative_groups.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <array>
#include <assert.h>
#include <atomic>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <stdint.h>
#include <vector>

#include "testdata.h"
#include "atomic.cuh"
#include "prefix.cuh"
#include "util.h"

namespace enumeration
{
__device__ __host__ inline unsigned int id() { return thread_id(); }

template <unsigned int levels> struct Buffer
{
  // shape [threads, levels]
  float *coefficients;
  // shape [threads, levels]
  float *partdists;
  // shape [threads, levels]
  float *centers;
  // shape [threads, levels]
  float *parent_center_partsums;

  constexpr static unsigned int size_in_bytes(unsigned int threads)
  {
    return 4 * sizeof(float) * threads * levels;
  }

  __device__ __host__ inline Buffer(unsigned char *memory, unsigned int threads)
      : coefficients(reinterpret_cast<float *>(memory)),
        partdists(reinterpret_cast<float *>(memory + sizeof(float) * threads * levels)),
        centers(reinterpret_cast<float *>(memory + 2 * sizeof(float) * threads * levels)),
        parent_center_partsums(
            reinterpret_cast<float *>(memory + 3 * sizeof(float) * threads * levels))
  {

  }

  __device__ __host__ inline float get_coeff(unsigned int level) const
  {
    return coefficients[id() * levels + level];
  }

  __device__ __host__ inline void set_coeff(unsigned int level, float coeff)
  {
    coefficients[id() * levels + level] = coeff;
  }

  __device__ __host__ inline float get_partdist(unsigned int level) const
  {
    return partdists[id() * levels + level];
  }

  __device__ __host__ inline void set_partdist(unsigned int level, float dist)
  {
    partdists[id() * levels + level] = dist;
  }

  __device__ __host__ inline float get_center(unsigned int level) const
  {
    return centers[id() * levels + level];
  }

  __device__ __host__ inline void set_center(unsigned int level, float dist)
  {
    centers[id() * levels + level] = dist;
  }

  __device__ __host__ inline float get_parent_center_partsum(unsigned int target) const
  {
    return parent_center_partsums[id() * levels + target];
  }

  __device__ __host__ inline void set_parent_center_partsum(unsigned int target, float dist)
  {
    parent_center_partsums[id() * levels + target] = dist;
  }
};

__device__ __host__ inline void next_coeff(float &coeff, const float center)
{
  const float rounded_center = round(center);
  if (isnan(coeff))
  {
    coeff = rounded_center;
  }
  else
  {
    coeff = 2 * rounded_center - coeff;
    if (center >= rounded_center)
    {
      coeff += static_cast<int>(coeff >= rounded_center);
    }
    else
    {
      coeff -= static_cast<int>(coeff <= rounded_center);
    }
  }
}

template <unsigned int levels, typename CallbackArg,
          void (*callback)(const Buffer<levels> &, float, unsigned int, CallbackArg&)>
__device__ __host__ inline void
diverging_enum(Buffer<levels> &buffer, const float *parent_points, const unsigned int parent_point_dimension,
               unsigned int *processed_parent_points, unsigned int parent_point_count, const float *mu, const unsigned int ldmu,
               const uint32_t *radius_square, unsigned long long *node_counter, CallbackArg& arg)
{
  unsigned int parent_index;
  unsigned int level = 0;

  float coefficient;
  float center;
  float partdist;
  while (true)
  {
    aggregated_atomic_inc(node_counter);
    const unsigned int k = levels - level - 1;
    if (level == 0)
    {
      parent_index = aggregated_atomic_inc(processed_parent_points);
      if (parent_index >= parent_point_count)
      {
        return;
      }
      for (unsigned int i = 0; i < levels; ++i)
      {
        float center_partsum = 0;
        for (unsigned int j = 0; j + 1 < parent_point_dimension; ++j)
        {
          center_partsum -=
              parent_points[parent_index * parent_point_dimension + j + 1] * mu[i * ldmu + j + levels];
        }
        buffer.set_parent_center_partsum(i, center_partsum);
      }
      partdist = 0;
      for (unsigned int j = 0; j + 1 < parent_point_dimension; ++j)
      {
        const float orthogonal_part = mu[(j + levels) * ldmu + j + levels];
        const float parent_coeff    = parent_points[parent_index * parent_point_dimension + j + 1];
        partdist += parent_coeff * parent_coeff * orthogonal_part *
                    orthogonal_part;
      }
      coefficient = parent_points[parent_index * parent_point_dimension];
      center      = buffer.get_parent_center_partsum(levels - 1) / mu[k * ldmu + k];
    }
    else
    {
      next_coeff(coefficient, center);
    }
    buffer.set_coeff(level, coefficient);
    buffer.set_center(level, center);
    buffer.set_partdist(level, partdist);

    const float alphak = coefficient - center;
    const float dist   = partdist + alphak * alphak * mu[k * ldmu + k] * mu[k * ldmu + k];
    const float limit  = int_to_float_order_preserving_bijection(*radius_square);

    assert(!isnan(limit));
    if (dist > limit)
    {
      if (level > 0)
      {
        --level;
      }
      coefficient = buffer.get_coeff(level);
      partdist    = buffer.get_partdist(level);
      center      = buffer.get_center(level);
    }
    else if (level == levels - 1)
    {
      callback(buffer, dist, parent_index, arg);
    }
    else
    {
      partdist = dist;
      center   = buffer.get_parent_center_partsum(k - 1);
      for (unsigned int i = 0; i <= level; ++i)
      {
        center -= buffer.get_coeff(i) * mu[(k - 1) * ldmu + levels - i - 1];
      }
      center /= mu[(k - 1) * ldmu + k - 1];
      coefficient = INFINITY - INFINITY;
      ++level;
    }
  }
}

struct CallbackData
{
  uint32_t *radius_squared_location;
  float *shortest_points_per_thread;
  const float *starting_points;
};

template <unsigned int levels>
__device__ __host__ inline void callback(const Buffer<levels> &buffer, float squared_norm, unsigned int parent_index,
                                         int &dummy)
{
  printf("Norm %f: ", squared_norm);
  float coefficient;
  for (unsigned int i = 0; i < levels; ++i)
  {
    coefficient = buffer.get_coeff(levels - i - 1);
    printf("%f, ", coefficient);
  }
  printf("\n");
}

template<unsigned int levels, unsigned int dimensions>
__device__ __host__ inline void process_point(const Buffer<levels> &buffer, float squared_norm,
                                              unsigned int parent_index, CallbackData &data)
{
  if (abs(squared_norm) < .5)
  {
    return;
  }
  float *shortest_points_per_thread = data.shortest_points_per_thread;
  uint32_t squared_norm_as_int = float_to_int_order_preserving_bijection(squared_norm);
  uint32_t old                 = atomic_min(data.radius_squared_location, squared_norm_as_int);
  if (squared_norm_as_int < old)
  {
    for (unsigned int i = 0; i < levels; ++i)
    {
      shortest_points_per_thread[thread_id() * dimensions + levels - i - 1] =
          buffer.get_coeff(i);
    }
    for (unsigned int i = 0; i < dimensions - levels; ++i)
    {
      shortest_points_per_thread[thread_id() * dimensions + levels + i] =
          data.starting_points[parent_index * (dimensions - levels + 1) + i + 1];
    }
  }
}

template <unsigned int levels, unsigned int dimensions>
__global__ void __launch_bounds__(512, 1) search_experiment(
      const float *starting_points, const unsigned int starting_point_count, unsigned int *processed_starting_points,
      unsigned char *buffer_memory, uint32_t *radius_squared_location, const float *mu, float *shortest_points_per_thread,
      unsigned long long *visited_node_counter)
{
  typedef cooperative_groups::thread_block CG;
  CG group = cooperative_groups::this_thread_block();

  CallbackData data = {radius_squared_location, shortest_points_per_thread, starting_points};
  Buffer<levels> buffer(buffer_memory, blockDim.x * gridDim.x);
  diverging_enum<levels, CallbackData, &process_point<levels, dimensions>>(
      buffer, starting_points, dimensions - levels + 1, processed_starting_points,
      starting_point_count, mu, dimensions, radius_squared_location, visited_node_counter, data);
}

template <unsigned int levels, unsigned int dimensions>
CudaPtr<float>
search_experiment_entry(const std::array<std::array<float, dimensions>, dimensions> &host_mu,
    const std::vector<std::array<float, dimensions - levels + 1>> &starting_points,
                   unsigned int &output_point_count, const float initial_radius)
{

  constexpr unsigned int block_size = 512;
  constexpr unsigned int grid_size  = 64;

  const uint32_t radius_store_init =
      float_to_int_order_preserving_bijection(initial_radius * initial_radius);
  CudaPtr<uint32_t> radius_squared = alloc(uint32_t, 1);
  check(cudaMemcpy(radius_squared.get(), &radius_store_init, sizeof(uint32_t),
                   cudaMemcpyHostToDevice));

  constexpr unsigned int starting_point_dim = dimensions - levels + 1;
  CudaPtr<float> device_starting_points = alloc(float, starting_points.size() * starting_point_dim);
  for (unsigned int i = 0; i < starting_points.size(); ++i)
  {
    check(cudaMemcpy(&device_starting_points.get()[i * starting_point_dim],
                     starting_points[i].data(), starting_point_dim * sizeof(float),
                     cudaMemcpyHostToDevice));
  }

  CudaPtr<unsigned int> processed_starting_points = alloc(unsigned int, 1);
  check(cudaMemset(processed_starting_points.get(), 0, sizeof(unsigned int)));

  CudaPtr<float> shortest_points_per_thread = alloc(float, grid_size * block_size * dimensions);

  CudaPtr<float> mu = alloc(float, dimensions * dimensions);
  for (unsigned int i = 0; i < dimensions; ++i)
  {
    check(cudaMemcpy(&mu.get()[i * dimensions], host_mu[i].data(), dimensions * sizeof(float),
                     cudaMemcpyHostToDevice));
  }

  std::cout << "Start " << grid_size << " blocks with " << block_size
            << " threads, global memory "
            << (Buffer<levels>::size_in_bytes(block_size * grid_size)) << " bytes"
            << std::endl;

  CudaPtr<unsigned char> buffer_memory =
      alloc(unsigned char, Buffer<levels>::size_in_bytes(block_size * grid_size));
  CudaPtr<unsigned long long> visited_node_counter = alloc(unsigned long long, 5);

  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

  search_experiment<levels, dimensions><<<dim3(grid_size), dim3(block_size)>>>(
              device_starting_points.get(), starting_points.size(), processed_starting_points.get(),
              buffer_memory.get(), radius_squared.get(), mu.get(), shortest_points_per_thread.get(),
              visited_node_counter.get());
  check(cudaDeviceSynchronize());
  check(cudaGetLastError());

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  uint32_t end_radius_square;
  check(cudaMemcpy(&end_radius_square, radius_squared.get(), sizeof(uint32_t),
                   cudaMemcpyDeviceToHost));
  std::cout << "End radius: " << sqrt(int_to_float_order_preserving_bijection(end_radius_square))
            << std::endl;

  unsigned long long counter;
  check(cudaMemcpy(&counter, visited_node_counter.get(), sizeof(unsigned long long),
                   cudaMemcpyDeviceToHost));
  std::cout << "Visited " << counter << " nodes in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
            << std::endl;

  print_performance_counter(&visited_node_counter.get()[1]);
  print_performance_counter(&visited_node_counter.get()[2]);
  print_performance_counter(&visited_node_counter.get()[3]);
  print_performance_counter(&visited_node_counter.get()[4]);

  output_point_count = grid_size * block_size;
  return shortest_points_per_thread;
}

inline void gpu_test() {
  constexpr unsigned int dimensions = 50;
  constexpr unsigned int levels     = dimensions - 2;

  const std::array<std::array<float, dimensions>, dimensions> &mu = test_mu_normal;

  std::vector<std::array<float, dimensions - levels + 1>> start_points;

  float initial_radius = find_initial_radius(mu) * 1.1;

  std::array<float, 3> x;
  std::function<void(const std::array<float, 3> &)> fn = [&start_points](auto point) {
    start_points.emplace_back(point);
  };
  CenterToOutIterator iter;
  do
  {
    x[2] = *iter;
    ++iter;
  } while (!naive_enum_recursive<3, dimensions>(x, 0.f, 0.f, mu, dimensions - 1,
                                                initial_radius * initial_radius, fn));

  unsigned int output_point_count           = 0;
  CudaPtr<float> shortest_points_per_thread = search_experiment_entry<levels, dimensions>(
      mu, start_points, output_point_count, find_initial_radius(mu) * 1.1);

  std::unique_ptr<float[]> result(new float[dimensions * output_point_count]);
  check(cudaMemcpy(result.get(), shortest_points_per_thread.get(),
                   dimensions * output_point_count * sizeof(float), cudaMemcpyDeviceToHost));

  if (!contains_solution(result.get(), test_solution_normal.data(), dimensions, output_point_count))
  {
    throw "Test failed";
  }
}

inline void cpu_test() {

  constexpr unsigned int dimensions = 20;
  constexpr unsigned int levels     = dimensions;

  std::vector<std::array<float, dimensions - levels + 1>> start_points;
  start_points.push_back({0});
  start_points.push_back({1});
  start_points.push_back({2});
  start_points.push_back({3});

  Buffer<levels> buffer(new unsigned char[Buffer<levels>::size_in_bytes(1)], 1);

  unsigned int processed_parent_points = 0;
  uint32_t radius_squared = float_to_int_order_preserving_bijection(
      find_initial_radius(test_mu_small) * 1.1 * find_initial_radius(test_mu_small) * 1.1);
  unsigned long long node_counter          = 0;
  int dummy;
  diverging_enum<levels, int, &callback>(buffer, &start_points[0][0], dimensions - levels + 1,
                                         &processed_parent_points, start_points.size(), & test_mu_small[0][0], dimensions,
                                         &radius_squared, &node_counter, dummy);
}

inline void test() { gpu_test(); }
}  // namespace enumeration
