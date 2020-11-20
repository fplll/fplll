#ifndef FPLLL_ATOMIC_CUH
#define FPLLL_ATOMIC_CUH

#include <assert.h>
#include <atomic>

#ifdef __CUDACC__
#define DEVICE_HOST_FUNCTION __device__ __host__
#include "cuda_runtime.h"
#include "cooperative_groups.h"
#else
#define DEVICE_HOST_FUNCTION
#endif

template<typename F, typename T> DEVICE_HOST_FUNCTION T bitcast(F bits)
{
    static_assert(sizeof(F) == sizeof(T), "Can only bit-cast types of same size");
    union Convert {
        T out;
        F in;
    };
    Convert result;
    result.in = bits;
    return result.out;
}

DEVICE_HOST_FUNCTION inline unsigned int float_to_int_order_preserving_bijection(float value)
{
    unsigned int val = bitcast<float, unsigned int>(value);
    unsigned int flip_all_if_negative_mask = static_cast<unsigned int>(-static_cast<int>(val >> 31));
    return val ^ (flip_all_if_negative_mask | 0x80000000);
}

DEVICE_HOST_FUNCTION inline float int_to_float_order_preserving_bijection(unsigned int value)
{
    unsigned int flip_all_if_negative_mask = static_cast<unsigned int>(-static_cast<int>((value >> 31) ^ 1));
    return bitcast<unsigned int, float>(value ^ (flip_all_if_negative_mask | 0x80000000));
}

// implementation is only atomic on the device
DEVICE_HOST_FUNCTION inline unsigned long long atomic_add(unsigned long long *target,
                                                         unsigned long long amount)
{
#ifdef __CUDA_ARCH__
    return atomicAdd(target, amount);
#else
    return (*target += amount) - amount;
#endif
}

// implementation is only atomic on the device, requires ctr to be the same for all threads in the warp
DEVICE_HOST_FUNCTION inline uint32_t aggregated_atomic_inc(uint32_t *ctr)
{
#ifdef __CUDA_ARCH__
  // use warp-aggregated atomic for improved performance
  auto g = cooperative_groups::coalesced_threads();
  uint32_t warp_res;
  if (g.thread_rank() == 0)
  {
    warp_res = atomicAdd(ctr, g.size());
  }
  return g.shfl(warp_res, 0) + g.thread_rank();
#else
  return (*ctr)++;
#endif
}

// implementation is only atomic on the device, requires ctr to be the same for all threads in the warp
DEVICE_HOST_FUNCTION inline uint64_t aggregated_atomic_inc(uint64_t *ctr)
{
  static_assert(sizeof(unsigned long long) == sizeof(uint64_t), "atomicAdd is only defined for unsigned long long, so we can only use it for uint64_t if unsigned long long == uint64_t");
#ifdef __CUDA_ARCH__
  // use warp-aggregated atomic for improved performance
  auto g = cooperative_groups::coalesced_threads();
  uint64_t warp_res;
  if (g.thread_rank() == 0)
  {
    warp_res = atomicAdd(reinterpret_cast<unsigned long long*>(ctr), static_cast<unsigned long long>(g.size()));
  }
  return g.shfl(warp_res, 0) + g.thread_rank();
#else
  return (*ctr)++;
#endif
}

// implementation is only atomic on the device
DEVICE_HOST_FUNCTION inline unsigned int atomic_add(unsigned int *target, unsigned int amount)
{
#ifdef __CUDA_ARCH__
  return atomicAdd(target, amount);
#else
  return (*target += amount) - amount;
#endif
}

// implementation is only atomic on the device
DEVICE_HOST_FUNCTION inline uint32_t atomic_min(uint32_t *target, uint32_t value)
{
#ifdef __CUDA_ARCH__
  return atomicMin(target, value);
#else
  uint32_t result = *target;
  if (value < *target)
  {
    *target = value;
  }
  return result;
#endif
}

// implementation is only atomic on the device
DEVICE_HOST_FUNCTION inline float atomic_min(uint32_t *target, float value)
{
  return int_to_float_order_preserving_bijection(
      atomic_min(target, float_to_int_order_preserving_bijection(value)));
}

// implementation is only atomic on the device
DEVICE_HOST_FUNCTION inline uint32_t atomic_load(volatile uint32_t *target)
{
  return *target;
}

DEVICE_HOST_FUNCTION inline void threadfence_system() {
#ifdef __CUDA_ARCH__
  __threadfence_system();
#endif
}

#endif