#ifndef FPLLL_PREFIX_CUH
#define FPLLL_PREFIX_CUH

#include "cuda_runtime.h"
#include "cooperative_groups.h"
#include <assert.h>

constexpr inline unsigned int int_log2(unsigned int x) {
    if (x == 0) {
        return std::numeric_limits<unsigned int>::max();
    }
    unsigned int result = 0;
    while (x >>= 1) {
        result += 1;
    }
    return result;
}

static_assert(int_log2(8) == 3, "error in log2");

template<typename CG, unsigned int block_size>
class PrefixCounter;

template<unsigned int block_size>
class PrefixCounter<cooperative_groups::thread_block, block_size> {

    unsigned char* shared_mem;

    /**
    Performs reduction from 2^(level + 1) sums to 2^level sums, and updates the cell_offsets accordingly.
    */
    template<unsigned int level>
    __device__ inline void prefix_sum_reduction_step(cooperative_groups::thread_block& group, unsigned int tid, unsigned int* accumulator, unsigned int* cell_offsets) {
        if constexpr (block_size_log >= level + 1) {
            unsigned int buffer;
            if (tid < (1 << level)) {
                buffer = accumulator[2 * tid] + accumulator[2 * tid + 1];
            }
            if ((tid >> (block_size_log - level - 1)) & 1) {
                cell_offsets[tid] += accumulator[(tid >> (block_size_log - level - 1)) - 1];
            }
            group.sync();
            if (tid < (1 << level)) {
                accumulator[tid] = buffer;
            }
            group.sync();
        }
    }

    template<unsigned int level>
    __device__ inline void prefix_count_reduction_step(cooperative_groups::thread_block& group, unsigned int tid, unsigned int* accumulator, unsigned int* cell_offsets) {
        constexpr unsigned int warpCountLog = block_size_log - 5;
        if constexpr (warpCountLog >= level + 1) {
            prefix_sum_reduction_step<level>(group, tid, accumulator, cell_offsets);
        }
    }

public:

    static_assert(block_size == (1L << int_log2(block_size)), "Expected BlockSize to be a power of 2");
    constexpr static unsigned int block_size_log = int_log2(block_size);

    __device__ inline PrefixCounter(unsigned char *shared_mem) : shared_mem(shared_mem) {}

    constexpr static unsigned int shared_mem_size_in_bytes = (block_size + block_size / 32) * sizeof(unsigned int);

    __device__ inline unsigned int prefix_count(cooperative_groups::thread_block &group,
                                                bool predicate, unsigned int &total_len)
    {
        assert(blockDim.x == block_size);
        assert(blockDim.y == 1 && blockDim.z == 1);

        unsigned int* cell_offsets = reinterpret_cast<unsigned int*>(shared_mem);
        unsigned int* accumulator = &reinterpret_cast<unsigned int*>(shared_mem)[block_size];

        const unsigned int warpid = threadIdx.x / 32;
        const unsigned int laneid = threadIdx.x % 32;
        const unsigned int in_warp_values = __ballot_sync(0xFFFFFFFF, predicate ? 1 : 0);
        const unsigned int in_warp_accumulation = __popc(in_warp_values);
        const unsigned int in_warp_prefix_count = __popc(in_warp_values << (32 - laneid));

        cell_offsets[threadIdx.x] = in_warp_prefix_count;
        if (laneid == 0) {
            accumulator[warpid] = in_warp_accumulation;
        }
        group.sync();

        const unsigned int tid = threadIdx.x;

        // now perform standard reduction with the remaining 1024/32 == 32 values
        prefix_count_reduction_step<4>(group, tid, accumulator, cell_offsets);
        prefix_count_reduction_step<3>(group, tid, accumulator, cell_offsets);
        prefix_count_reduction_step<2>(group, tid, accumulator, cell_offsets);
        prefix_count_reduction_step<1>(group, tid, accumulator, cell_offsets);
        prefix_count_reduction_step<0>(group, tid, accumulator, cell_offsets);

        group.sync();

        total_len = accumulator[0];
        return cell_offsets[threadIdx.x];
    }
};

template <unsigned int block_size> class PrefixCounter<cooperative_groups::thread_block_tile<32>, block_size>
{

public:
  static_assert(block_size == (1L << int_log2(block_size)),
                "Expected BlockSize to be a power of 2");
  constexpr static unsigned int block_size_log = int_log2(block_size);

  __device__ inline PrefixCounter() {}

  constexpr static unsigned int shared_mem_size_in_bytes = 0;

  __device__ inline unsigned int prefix_count(cooperative_groups::thread_block_tile<32> &group,
                                              bool predicate, unsigned int &total_len)
  {
    assert(blockDim.x == block_size);
    assert(blockDim.y == 1 && blockDim.z == 1);

    const unsigned int warpid               = threadIdx.x / 32;
    const unsigned int laneid               = threadIdx.x % 32;
    const unsigned int in_warp_values       = __ballot_sync(0xFFFFFFFF, predicate ? 1 : 0);
    const unsigned int in_warp_accumulation = __popc(in_warp_values);
    const unsigned int in_warp_prefix_count = __popc(in_warp_values << (32 - laneid));

    total_len = in_warp_accumulation;
    return in_warp_prefix_count;
  }
};

struct single_thread {

    __device__ __host__ inline void sync() {}

    __device__ __host__ inline unsigned int thread_rank() {
        return 0;
    }

    __device__ __host__ inline unsigned int size() {
        return 1;
    }
};

template<unsigned int block_size>
class PrefixCounter<single_thread, block_size> {

public:

    static_assert(block_size == 1, "Using single_thread group requires the block size to be exactly 1 thread");

    __device__ __host__ inline unsigned int prefix_count(single_thread &group, bool predicate,
                                                         unsigned int &total_len)
    {
        total_len = predicate ? 1 : 0;
        return 0;
    }
};

#endif