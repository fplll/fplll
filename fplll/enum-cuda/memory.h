#ifndef FPLLL_CUDA_MEMORY_H
#define FPLLL_CUDA_MEMORY_H

#include "cuda_runtime.h"
#include <memory>
#include "cuda_check.h"

template <typename T> struct CudaDeleter
{
  inline void operator()(T *ptr) const { check(cudaFree(ptr)); }
};

template <typename T> struct PinnedDeleter
{
  inline void operator()(T *ptr) const { check(cudaFreeHost(ptr)); }
};

struct EventDeleter
{
  inline void operator()(cudaEvent_t ptr) const { check(cudaEventDestroy(ptr)); }
};

struct StreamDeleter
{
  inline void operator()(cudaStream_t ptr) const { check(cudaStreamDestroy(ptr)); }
};

template<typename T>
using CudaPtr = std::unique_ptr<T, CudaDeleter<T>>;

template<typename T>
using PinnedPtr = std::unique_ptr<T, PinnedDeleter<T>>;

using CudaEvent = std::unique_ptr<std::remove_pointer<cudaEvent_t>::type, EventDeleter>;

using CudaStream = std::unique_ptr<std::remove_pointer<cudaStream_t>::type, StreamDeleter>;

template<typename T>
inline CudaPtr<T> allocateCudaMemory(size_t len, const char* file, const unsigned int line)
{
	T* result = nullptr;
	checkCudaError(cudaMalloc(&result, len * sizeof(T)), file, line);
	return CudaPtr<T>(result);
}

template<typename T>
inline PinnedPtr<T> allocatePinnedMemory(size_t len)
{
	T* result = nullptr;
	checkCudaError(cudaMallocHost(&result, len * sizeof(T)), __FILE__, __LINE__);
	return PinnedPtr<T>(result);
}

#define cuda_alloc(type, len) allocateCudaMemory<type>(len, __FILE__, __LINE__) 



#endif