#pragma once
#include "cuda_runtime.h"
#include "cooperative_groups.h"
#include <assert.h>
#include <atomic>

template<typename F, typename T>
__device__ __host__ T bitcast(F bits) {
    static_assert(sizeof(F) == sizeof(T), "Can only bit-cast types of same size");
    union Convert {
        T out;
        F in;
    };
    Convert result;
    result.in = bits;
    return result.out;
}

__device__ __host__ inline unsigned int float_to_int_order_preserving_bijection(float value)
{
    unsigned int val = bitcast<float, unsigned int>(value);
    unsigned int flip_all_if_negative_mask = static_cast<unsigned int>(-static_cast<int>(val >> 31));
    return val ^ (flip_all_if_negative_mask | 0x80000000);
}

__device__ __host__ inline float int_to_float_order_preserving_bijection(unsigned int value)
{
    unsigned int flip_all_if_negative_mask = static_cast<unsigned int>(-static_cast<int>((value >> 31) ^ 1));
    return bitcast<unsigned int, float>(value ^ (flip_all_if_negative_mask | 0x80000000));
}

__device__ __host__ inline unsigned long long atomic_add(unsigned long long *target,
                                                         unsigned long long amount)
{
#ifdef __CUDA_ARCH__
    return atomicAdd(target, amount);
#else
    return (*target += amount) - amount;
#endif
}

__device__ __host__ inline uint32_t atomic_min(uint32_t* target, uint32_t value) {
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

__device__ __host__ inline float atomic_min(uint32_t *target, float value) {
  return int_to_float_order_preserving_bijection(
      atomic_min(target, float_to_int_order_preserving_bijection(value)));
}