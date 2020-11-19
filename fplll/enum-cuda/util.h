#pragma once
#include <iostream>
#include <functional>
#include <cmath>
#include <vector>

/**
Provides an iterator that iterates over Z in the order 0, -1, 1, -2, 2, -3, 3, ...
*/
struct CenterToOutIterator
{
  int current;

  inline CenterToOutIterator() : current(0) {}

  inline void operator++()
  {
    current += static_cast<int>(current >= 0);
    current = -current;
  }

  inline void operator+=(unsigned int val)
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

  inline int operator*() { return current; }
};

template <typename FT, typename FL>
inline bool naive_enum_recursive(std::vector<FT> &x, 
                                 unsigned int levels, unsigned int n,
                                 const FL parentdist, const FL parentcenter,
                                 const FL* mu,
                                 const unsigned int k, const FL radius_square,
                                 std::function<void(const std::vector<FT> &)> &callback)
{
  FL alphak  = x[k + levels - n] - parentcenter;
  FL newdist = parentdist + alphak * alphak * mu[k * n + k] * mu[k * n + k];

  if (newdist > radius_square)
  {
    return true;
  }

  if (k == n - levels)
  {
    callback(x);
    return false;
  }

  FL newcenter = 0;
  for (unsigned int i = k; i < levels; ++i)
  {
    newcenter -= x[i] * mu[(k - 1) * n + i];
  }
  newcenter /= mu[(k - 1) * n + k - 1];

  bool is_out_of_bounds = false;
  for (CenterToOutIterator iter; !is_out_of_bounds; ++iter)
  {
    x[k + levels - n - 1] = round(newcenter) + *iter;
    is_out_of_bounds =
        naive_enum_recursive<FT>(x, levels, n, newdist, newcenter, mu, k - 1, radius_square, callback);
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
    result = std::min<FL>(result, norm_square);
  }
  return sqrt(result);
}

struct FloatWrapper {
  float value;

  inline double get_d() const {
    return value;
  }

  inline FloatWrapper operator-(float rhs) const {
    return { value - rhs };
  }

  inline FloatWrapper operator*(float rhs) const {
    return { value * rhs };
  }

  inline FloatWrapper& operator=(float rhs) {
    value = rhs;
    return *this;
  }

  inline FloatWrapper& operator++() {
    ++value;
    return *this;
  }

  inline operator float() const {
    return value;
  }
};
