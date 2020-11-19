#ifndef FPLLL_STREAMING_CUH
#define FPLLL_STREAMING_CUH

#include "cooperative_groups.h"
#include "cuda_runtime.h"
#include <assert.h>
#include <cstring>
#include <vector>
#include <iostream>

#include "atomic.h"
#include "memory.h"
#include "cuda_check.h"
#include "cuda_util.cuh"
#include "types.cuh"

/**
 * Evaluator proxy that streams all found points to the host, via a circular buffer
 */
template<unsigned int buffer_size>
class PointStreamEvaluator
{
  // these are accessible by the host
  enumi* result_buffer;
  enumf* result_squared_norms;
  unsigned int* has_written_round;
  
  // these are internal for collecting the points
  unsigned int* write_to_index;
  unsigned int point_index;
  
public:
  
  __device__ __host__ inline PointStreamEvaluator(unsigned char* memory, unsigned int* write_to_index)
    : result_squared_norms(reinterpret_cast<enumf*>(memory)),
      has_written_round(reinterpret_cast<unsigned int*>(memory + sizeof(enumf) * buffer_size)),
      result_buffer(reinterpret_cast<enumi*>(memory + (sizeof(enumf) + sizeof(unsigned int)) * buffer_size)),
      write_to_index(write_to_index)
  {
  }
 
  __device__ __host__ constexpr static unsigned int memory_size_in_bytes(unsigned int point_dimension) {
    unsigned int real_size = sizeof(enumi) * buffer_size * point_dimension + (sizeof(enumf) + sizeof(unsigned int)) * buffer_size;
    // ensure alignment
    return ((real_size - 1) / sizeof(enumf) + 1) * sizeof(enumf);
  }

  __device__ __host__ inline void operator()(enumi x, unsigned int coordinate, bool done, enumf norm_square,
                                             uint32_t *enum_bound_location)
  {
    if (coordinate == 0) {
      point_index = aggregated_atomic_inc(write_to_index) % buffer_size;
      result_squared_norms[point_index] = norm_square;
    }
 
    result_buffer[point_index + coordinate * buffer_size] = x;
 
    if (done) {
      threadfence_system();
      atomic_add(&has_written_round[point_index], 1);
    }
  }
};

template<unsigned int buffer_size>
class PointStreamEndpoint {

  void* device_memory;
  uint32_t* device_enumeration_bound_location;
  unsigned int evaluator_count;
  unsigned int point_dimension;
  PinnedPtr<unsigned char> host_memory;
  std::vector<std::unique_ptr<unsigned int[]>> last_has_written_round;
  CudaStream stream;

  inline void* evaluator_memory(unsigned int evaluator_id) {
      return static_cast<void*>(host_memory.get() + evaluator_id * PointStreamEvaluator<buffer_size>::memory_size_in_bytes(point_dimension));
  }

  inline unsigned int* host_has_written_round(unsigned int evaluator_id) {
    return static_cast<unsigned int*>(evaluator_memory(evaluator_id) + sizeof(enumf) * buffer_size);
  }

  inline unsigned int* device_has_written_round(unsigned int evaluator_id) {
    void* device_evaluator_memory = device_memory + PointStreamEvaluator<buffer_size>::memory_size_in_bytes(point_dimension) * evaluator_id;
    return static_cast<unsigned int*>(device_evaluator_memory + sizeof(enumf) * buffer_size);
  }
  
public:

  inline PointStreamEndpoint(unsigned char* device_memory, uint32_t* device_enumeration_bound_location, unsigned int evaluator_count, unsigned int point_dimension)
    : evaluator_count(evaluator_count),
      point_dimension(point_dimension),
      device_enumeration_bound_location(device_enumeration_bound_location),
      device_memory(device_memory),
      host_memory(allocatePinnedMemory<unsigned char>(evaluator_count * PointStreamEvaluator<buffer_size>::memory_size_in_bytes(point_dimension)))
    {
        const unsigned int used_device = 0;
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, used_device);
        if (deviceProp.asyncEngineCount <= 1) {
            throw "your cuda device does not support asynchronous kernel execution and data transfer, which is required for solution point streaming";
        }
        cudaStream_t raw_stream;
        check(cudaStreamCreate(&raw_stream));
        stream.reset(raw_stream);
        for (unsigned int i = 0; i < evaluator_count; ++i) {
            last_has_written_round.push_back(std::unique_ptr<unsigned int[]>(new unsigned int[buffer_size]));
        }
    }

    __host__ inline void init() {
        for (unsigned int i = 0; i < evaluator_count; ++i) {
            // use default cuda stream to prevent overlap with kernel execution & point query
            check(cudaMemset(device_has_written_round(i), 0, buffer_size * sizeof(unsigned int)));
            std::memset(host_has_written_round(i), 0, buffer_size * sizeof(unsigned int));
        }
    }
  
    template<typename Fn, bool print_status = true>
    __host__ inline void query_new_points(Fn callback) {
        for (unsigned int i = 0; i < evaluator_count; ++i) {
            std::memcpy(last_has_written_round[i].get(), host_has_written_round(i), buffer_size * sizeof(unsigned int));
        }
        const unsigned int total_size = PointStreamEvaluator<buffer_size>::memory_size_in_bytes(point_dimension) * evaluator_count;
        check(cudaMemcpyAsync(host_memory.get(), device_memory, total_size, cudaMemcpyDeviceToHost, stream.get()));
        check(cudaStreamSynchronize(stream.get()));

        float new_enumeration_bound = INFINITY;
        for (unsigned int i = 0; i < evaluator_count; ++i) {
            for (unsigned int j = 0; j < buffer_size; ++j) {
                const unsigned int last_written_count = last_has_written_round[i][j];
                const unsigned int new_written_count = host_has_written_round(i)[j];
                if (last_written_count == new_written_count) {
                    // no new data here in this buffer entry
                } else if (last_written_count + 1 == new_written_count) {
                    enumf norm_square = static_cast<enumf*>(evaluator_memory(i))[j];
                    enumi* x = static_cast<enumi*>(evaluator_memory(i) + (sizeof(enumf) + sizeof(unsigned int)) * buffer_size);
                    enumf new_enum_bound = callback(norm_square, x);
                    new_enumeration_bound = std::min<float>(new_enum_bound, new_enumeration_bound);
                } else {
                    throw "buffer not big enough to hold all solution points found between two queries";
                }
            }
        }
        if (!isinf(new_enumeration_bound)) {
            uint32_t enumeration_bound_encoded = float_to_int_order_preserving_bijection(new_enumeration_bound);
            check(cudaMemcpyAsync(device_enumeration_bound_location, &enumeration_bound_encoded, sizeof(uint32_t), cudaMemcpyHostToDevice, stream.get()));
            if (print_status) {
                std::cout << "Can decrease enum bound to " << new_enumeration_bound << std::endl;
            }
        }
    } 

    __host__ inline void wait_for_event(cudaEvent_t event) {
        check(cudaStreamWaitEvent(stream.get(), event, 0));
    }
};

#endif