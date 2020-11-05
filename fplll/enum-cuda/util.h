#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <functional>

/**
Provides an iterator that iterates over Z in the order 0, -1, 1, -2, 2, -3, 3, ...
*/
struct CenterToOutIterator
{
  int current;

  __device__ __host__ inline CenterToOutIterator() : current(0) {}

  __device__ __host__ inline void operator++()
  {
    current += static_cast<int>(current >= 0);
    current = -current;
  }

  __device__ __host__ inline void operator+=(unsigned int val)
  {
    if (val % 2 == 0)
    {
      current += current >= 0 ? static_cast<int>(val) / 2 : -static_cast<int>(val) / 2;
    }
    else
    {
      operator+=(val - 1);
      operator++();
    }
  }

  __device__ __host__ inline int operator*() { return current; }
};

template <unsigned int levels, unsigned int dimensions>
inline bool naive_enum_recursive(std::array<float, levels> &x, const float parentdist,
                                 const float parentcenter,
                                 const std::array<std::array<float, dimensions>, dimensions> &mu,
                                 const unsigned int k, const float radius_square,
                                 std::function<void(const std::array<float, levels> &)> &callback)
{
  float alphak  = x[k + levels - dimensions] - parentcenter;
  float newdist = parentdist + alphak * alphak * mu[k][k] * mu[k][k];

  if (newdist > radius_square)
  {
    return true;
  }

  if (k == dimensions - levels)
  {
    callback(x);
    return false;
  }

  float newcenter = 0;
  for (unsigned int i = k; i < levels; ++i)
  {
    newcenter -= x[i] * mu[k - 1][i];
  }
  newcenter /= mu[k - 1][k - 1];

  bool is_out_of_bounds = false;
  for (CenterToOutIterator iter; !is_out_of_bounds; ++iter)
  {
    x[k + levels - dimensions - 1] = round(newcenter) + *iter;
    is_out_of_bounds =
        naive_enum_recursive(x, newdist, newcenter, mu, k - 1, radius_square, callback);
  }
  return false;
}

inline bool contains_solution(const float *shortest_vectors, const float *expected,
                       const unsigned int levels, const unsigned int vector_count)
{
  for (unsigned int vector_id = 0; vector_id < vector_count; ++vector_id)
  {
    bool contains_solution     = true;
    bool contains_neg_solution = true;
    for (unsigned int i = 0; i < levels; ++i)
    {
      // equality is ok, as shortest_vectors contains integers
      contains_solution &= shortest_vectors[vector_id * levels + i] == expected[i];
      contains_neg_solution &= shortest_vectors[vector_id * levels + i] == -expected[i];
    }
    if (contains_solution || contains_neg_solution)
    {
      return true;
    }
  }
  std::cout << "Expected" << std::endl;
  for (unsigned int i = 0; i < levels; ++i)
  {
    std::cout << expected[i] << ", ";
  }
  std::cout << std::endl << "Actual" << std::endl;
  for (unsigned int vector_id = 0; vector_id < vector_count; ++vector_id)
  {
    for (unsigned int i = 0; i < levels; ++i)
    {
      std::cout << shortest_vectors[vector_id * levels + i] << ", ";
    }
    std::cout << std::endl;
  }
  return false;
}

template <typename FL, unsigned int dimensions>
FL find_initial_radius(const std::array<std::array<FL, dimensions>, dimensions> &mu)
{
  FL result = INFINITY;
  for (unsigned int vector = 0; vector < dimensions; ++vector)
  {
    FL norm_square = 0;
    for (unsigned int i = 0; i <= vector; ++i)
    {
      norm_square += mu[i][vector] * mu[i][vector];
    }
    result = min(result, norm_square);
  }
  return sqrt(result);
}

inline void print_performance_counter(unsigned long long *device_ptr)
{
  unsigned long long counter;
  check(cudaMemcpy(&counter, device_ptr, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
  std::cout << "Performance counter: " << counter << std::endl;
}

template <typename CG>
__device__ __host__ bool all_threads_eq(CG &group, unsigned int value, unsigned int *shared_mem)
{
  group.sync();
  if (group.thread_rank() == 0)
  {
    *shared_mem = value;
  }
  group.sync();
  return value == *shared_mem;
}