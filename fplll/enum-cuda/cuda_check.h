#ifndef FPLLL_CUDA_CHECK_H
#define FPLLL_CUDA_CHECK_H

#include "cuda_runtime.h"
#include <iostream>

struct CudaError
{
  cudaError_t status;

  constexpr inline CudaError(cudaError_t status) : status(status) {}
};

inline void checkCudaError(cudaError_t status, const char *file, const unsigned int line)
{
  if (status != cudaSuccess)
  {
    std::cerr << "Error: " << cudaGetErrorString(status) << " at " << file << ":" << line
              << std::endl;
    throw CudaError(status);
  }
}

#define check(x) checkCudaError(x, __FILE__, __LINE__)

#endif
