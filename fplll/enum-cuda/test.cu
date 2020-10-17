#include "enum.cu"
#include "testdata.h"

bool contains_solution(const float* shortest_vectors, const float* expected, const unsigned int levels, const unsigned int vector_count) {
    for (unsigned int vector_id = 0; vector_id < vector_count; ++vector_id) 
    {
        bool contains_solution = true;
        bool contains_neg_solution = true;
        for (unsigned int i = 0; i < levels; ++i) 
        {
            // equality is ok, as shortest_vectors contains integers
            contains_solution &= shortest_vectors[vector_id * levels + i] == expected[i];
            contains_neg_solution &= shortest_vectors[vector_id * levels + i] == -expected[i];
        }
        if (contains_solution || contains_neg_solution) {
            return true;
        }
    }
    std::cout << "Expected" << std::endl;
    for (unsigned int i = 0; i < levels; ++i)
    {
      std::cout << expected[i] << ", ";
    }
    std::cout << std::endl << "Actual" << std::endl;
    for (unsigned int vector_id = 0; vector_id < vector_count; ++vector_id)
    {
      for (unsigned int i = 0; i < levels; ++i)
      {
        std::cout << shortest_vectors[vector_id * levels + i] << ", ";
      }
      std::cout << std::endl;
    }
    return false;
}

template<typename FL>
FL find_initial_radius(const FL* mu, unsigned int levels) {
    FL result = INFINITY;
    for (unsigned int vector = 0; vector < levels; ++vector) {
        FL norm_square = 0;
        for (unsigned int i = 0; i <= vector; ++i) {
            norm_square += mu[i + vector * levels] * mu[i + vector * levels];
        }
        result = min(result, norm_square);
    }
    return sqrt(result);
}

void simple_gpu_test() {

    constexpr unsigned int levels = 3;
    std::array<float, levels * levels> host_mu = { 4.58257569, 0., 0., 4.3643578, 1.71824939, 0., 3.27326835, 2.16166858, 1.27000127 };

    unsigned int output_point_count = 0;
    CudaPtr<float> shortest_points_per_thread = search_enumeration<float, levels>(host_mu.data(), output_point_count, find_initial_radius(host_mu.data(), levels) * 1.1);

    std::unique_ptr<float[]> result(new float[levels * output_point_count]);
    check(cudaMemcpy(result.get(), shortest_points_per_thread.get(), levels * output_point_count * sizeof(float), cudaMemcpyDeviceToHost));

    std::array<float, levels> expected = { 1, -1, 0 };
    if (!contains_solution(result.get(), expected.data(), levels, output_point_count)) {
        throw "Test failed";
    }
}

void complex_gpu_test() {

    constexpr unsigned int levels = 50;

    std::array<float, levels * levels> host_mu;
    for (unsigned int i = 0; i < levels; ++i) {
        for (unsigned int j = 0; j < levels; ++j) {
            host_mu[i + j * levels] = test_mu_normal[i][j];
        }
    }

    unsigned int output_point_count = 0;
    CudaPtr<float> shortest_points_per_thread =
        search_enumeration<float, levels>(host_mu.data(), output_point_count, find_initial_radius(host_mu.data(), levels) * 1.1);

    std::unique_ptr<float[]> result(new float[levels * output_point_count]);
    check(cudaMemcpy(result.get(), shortest_points_per_thread.get(), levels * output_point_count * sizeof(float), cudaMemcpyDeviceToHost));

    if (!contains_solution(result.get(), test_solution_normal.data(), levels, output_point_count)) {
        throw "Test failed";
    }
}

int main()
{
  simple_gpu_test();
  complex_gpu_test();

  return 0;
}
