#ifndef FPLLL_GPUFUNC_CUH
#define FPLLL_GPUFUNC_CUH

#include <memory>
#include <utility>

template<typename F>
class opaque_unified_function {
	std::unique_ptr<void, void(*)(void*)> the_function;

	inline opaque_unified_function(std::unique_ptr<void, void(*)(void*)>&& the_function) : the_function(std::move(the_function)) {}
};

#ifdef __CUDACC__

#include "cuda_runtime.h"
#include <type_traits>
#include <cstdlib>
#include <functional>
#include "cuda_util.cuh"

template<typename F>
class gpu_function;

template<typename F>
class unified_function;

/**
 * Stores a function that can be called from device code. The function does not
 * have to be known at compile time and may have (mutable) internal state. gpu_function
 * does not own its internal state.
 */
template<typename R, typename... Args>
class gpu_function<R(Args...)> {
	void* state;
	R (*function_ptr)(void*, Args...);

	gpu_function(void* state, R(*function_ptr)(void*, Args...)) : state(state), function_ptr(function_ptr) {}

public:

	__device__ inline void operator()(Args... args) {
		function_ptr(state, args...);
	}

	template<typename F>
	friend class unified_function;
};

/**
 * Stores a function that can be called both from host and from device code. The function does not
 * have to be known at compile time and may have (mutable) internal state. This state is not owned
 * by unified_function.
 */
template<typename R, typename... Args>
class unified_function<R(Args...)> {

	typedef R(*fn_t)(void*, Args...);

	void* state;
	void* gpu_state;
	unsigned int state_size_in_bytes;
	fn_t device_function_ptr;
	fn_t host_function_ptr;

	/** 
	 * Getting a correct device function pointer is nontrivial, therefore this constructor is private. Instead,
	 * the creator functions (e.g. unified_function::bind()) should be used instead.
	 */
	inline unified_function(void* state, unsigned int state_size_in_bytes, fn_t device_function_ptr, fn_t host_function_ptr)
		: state(state), state_size_in_bytes(state_size_in_bytes), 
		  device_function_ptr(device_function_ptr), 
		  host_function_ptr(host_function_ptr) 
	{
		check(cudaMalloc(&gpu_state, state_size_in_bytes));
	}

	static inline void delete_as_void_ptr(void* ptr) {
		delete static_cast<unified_function<R(Args...)>*>(ptr);
	}

public:

	inline unified_function(opaque_unified_function<R(Args...)>&& that) {
		unified_function<R(Args...)>* ptr = static_cast<unified_function<R(Args...)>*>(that.the_function.get());
		*this = std::move(*ptr);
	}

	inline operator opaque_unified_function<R(Args...)>() && {
		unified_function<R(Args...)>* ptr = new unified_function<R(Args...)>(std::move(*this));
		std::unique_ptr<void, void(*)(void*)> the_function(ptr, delete_as_void_ptr);
		return opaque_unified_function<R(Args...)>(std::move(the_function));
	}

	explicit inline unified_function(const unified_function& that) 
		: unified_function(that.state, that.state_size_in_bytes, that.device_function_ptr, that.host_function_ptr)
	{

	}

	inline unified_function(unified_function&& that)
	{
		*this = std::move(that);
	}

	inline ~unified_function() {
		check(cudaFree(gpu_state));
	}

	unified_function& operator=(const unified_function&) = delete;
	inline unified_function& operator=(unified_function&& that) {
		state = that.state;
		gpu_state = that.gpu_state;
		state_size_in_bytes = that.state_size_in_bytes;
		device_function_ptr = that.device_function_ptr;
		host_function_ptr = that.host_function_ptr;
		that.gpu_state = nullptr;
	}

	inline gpu_function<R(Args...)> get_device_function() {
		assert(gpu_state != nullptr);
		assert(device_function_ptr != nullptr);
		return gpu_function<R(Args...)>(gpu_state, device_function_ptr);
	}

	inline std::function<R(Args...)> get_host_function() {
		return [state, host_function_ptr](Args... args) -> R { return host_function_ptr(state, args...); };
	}

	inline void copy_state_to_device() {
		assert(gpu_state != nullptr);
		check(cudaMemcpy(gpu_state, state, state_size_in_bytes, cudaMemcpyHostToDevice));
	}

	inline void copy_state_to_host() {
		assert(gpu_state != nullptr);
		check(cudaMemcpy(state, gpu_state, state_size_in_bytes, cudaMemcpyDeviceToHost));
	}

	inline void operator()(Args... args) {
		host_function_ptr(state, args...);
	}

	template<typename T, R(T::* f)(Args...)>
	static unified_function bind(T& value);
};

namespace gpufunc_internal {

	template<typename T, typename R, typename... Args>
	struct MemberFunctionWrapper {
		typedef R(T::* member_function)(Args...);
		typedef R(*raw_function)(void*, Args...);
		typedef T class_type;

		template<member_function f>
		__device__ __host__ static R call_converted_member(void* state, Args... args) {
			return (static_cast<T*>(state)->*f)(args...);
		}
	};

	template<typename F, F::member_function f>
	__device__ F::raw_function generic_device_function_ptr = &F::template call_converted_member<f>;

	template<typename F, F::member_function f>
	constexpr F::raw_function generic_host_function_ptr = &F::template call_converted_member<f>;

}

/**
 * Creates a unified function corresponding to a __host__ __device__ member function on an object. Calling the resulting function will 
 * call the member function on the given object, therefore the reference to the object must be valid for the whole lifetime of the
 * returned unified_function. Using the function on the gpu requires copying the object state to gpu memory via memcpy, therefore the object
 * must be trivially copyable.
 */
template<typename R, typename... Args>
template<typename T, R(T::*f)(Args...)>
inline unified_function<R(Args...)> unified_function<R(Args...)>::bind(T& value) {
	static_assert(std::is_trivially_copyable_v<T>, "State of a gpu function must be mem-copyable to and from device memory");

	fn_t device_fn;
	check(cudaMemcpyFromSymbol(&device_fn, gpufunc_internal::generic_device_function_ptr<gpufunc_internal::MemberFunctionWrapper<T, R, Args...>, f>, sizeof(fn_t)));

	return unified_function(static_cast<void*>(&value), sizeof(T), device_fn, gpufunc_internal::generic_host_function_ptr<gpufunc_internal::MemberFunctionWrapper<T, R, Args...>, f>);
}

template<typename F>
using unified_function_ptr = std::shared_ptr<unified_function<F>>;

#endif

#endif