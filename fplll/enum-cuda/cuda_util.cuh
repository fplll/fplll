#ifndef FPLLL_CUDA_UTIL_CUH
#define FPLLL_CUDA_UTIL_CUH

#include "cuda_runtime.h"
#include <iostream>
#include <memory>
#include <ctime>
#include "memory.h"

__device__ __host__ inline unsigned int thread_id()
{
#ifdef __CUDA_ARCH__
  return threadIdx.x + blockIdx.x * blockDim.x;
#else
  return 0;
#endif
}

__device__ __host__ inline unsigned int thread_id_in_block()
{
#ifdef __CUDA_ARCH__
  return threadIdx.x;
#else
  return 0;
#endif
}

__device__ __host__ inline unsigned long long time()
{
#ifdef __CUDA_ARCH__
  return clock64();
#else
  return clock();
#endif
}

__device__ __host__ inline void runtime_error() {
#ifdef __CUDA_ARCH__
  asm("trap;");
#else
  throw;
#endif
}

#endif