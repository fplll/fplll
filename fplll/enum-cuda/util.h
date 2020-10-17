#pragma once
#include "cuda_runtime.h"
#include <iostream>

constexpr bool TRACE = false;

struct CudaError {
	cudaError_t status;

	constexpr inline CudaError(cudaError_t status) : status(status)
	{
	}
};

inline void checkCudaError(cudaError_t status, const char *file, const unsigned int line)
{
	if (status != cudaSuccess) {
		std::cerr << "Error: " << cudaGetErrorString(status) << " at " << file << ":" << line << std::endl;
		throw CudaError(status);
	}
}

#define check(x) checkCudaError(x, __FILE__, __LINE__)

template<typename T>
class CudaDeleter
{
public:
	CudaDeleter() noexcept = default;
	CudaDeleter(const CudaDeleter&) noexcept = default;
	CudaDeleter(CudaDeleter&&) noexcept = default;
	~CudaDeleter() = default;

	CudaDeleter& operator=(const CudaDeleter&) = default;
	CudaDeleter& operator=(CudaDeleter&&) = default;

	void operator()(T*) const;
};

template<typename T>
inline void CudaDeleter<T>::operator()(T* ptr) const
{
	check(cudaFree(ptr));
}

template<typename T>
using CudaPtr = std::unique_ptr<T, CudaDeleter<T>>;

template<typename T>
inline CudaPtr<T> allocateCudaMemory(size_t len, const char* file, const unsigned int line)
{
	T* result = nullptr;
	checkCudaError(cudaMalloc(&result, len * sizeof(T)), file, line);
	return CudaPtr<T>(result);
}

#define alloc(type, len) allocateCudaMemory<type>(len, __FILE__, __LINE__) 

__device__ __host__ inline unsigned int thread_id()
{
#ifdef __CUDA_ARCH__
  return threadIdx.x + blockIdx.x * blockDim.x;
#else
  return 0;
#endif
}

template <typename... Args> void debug_message_thread(const char *string, Args... arguments)
{
  if (thread_id() == 0 && TRACE)
  {
    printf(string, arguments...);
  }
}