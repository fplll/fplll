#ifndef FPLLL_TYPES_CUH
#define FPLLL_TYPES_CUH

#include "cuda_runtime.h"
#include "atomic.h"

typedef double enumf;
typedef float enumi;

struct Matrix
{
  const enumf *ptr;
  unsigned int ld;

  __device__ __host__ inline Matrix(const enumf *ptr, unsigned int ld) : ptr(ptr), ld(ld) {}

  __device__ __host__ inline Matrix() : ptr(nullptr), ld(0) {}

  __device__ __host__ inline enumf at(unsigned int row, unsigned int col) const
  {
    return ptr[row * ld + col];
  }

  __device__ __host__ inline Matrix block(unsigned int start_row, unsigned int start_col) const
  {
    return Matrix(&ptr[start_col + start_row * ld], ld);
  }
};

template <unsigned int levels, unsigned int dimensions_per_level, unsigned int max_nodes_per_level>
struct SubtreeEnumerationBuffer;

struct PerfCounter
{
  unsigned long long *counter;

  __device__ __host__ inline PerfCounter(unsigned long long *target) : counter(target) {}

  __device__ __host__ inline void inc() { atomic_inc(counter); }
};

#endif