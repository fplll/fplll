#ifndef FPLLL_CUDA_ENUM
#define FPLLL_CUDA_ENUM

enum struct EvaluatorStrategy
{
    FAST
};

struct CudaEnumOpts
{
  // maximal amount of paths that will be searched using recursive enumeration during each algorithm step
  // by each thread. When this is exceeded, the recursive enumeration state is stored and resumed in the 
  // next step. Use for load balancing to prevent threads with small subtrees to wait for threads with very big subtrees.
  unsigned int max_subtree_paths;
  // stop children generation when the percentage of parent points that have still unprocessed
  // children drops beneath this percentage
  float min_active_parents_percentage;
  // height of the subtrees that are searched using recursive enumeration.
  unsigned int dimensions_per_level;
  // how many start points should be assigned to each cooperative group
  unsigned int initial_nodes_per_group;
  // how many cuda threads to use for the search. If this is not a multiple of the block
  // size, it will be rounded up to one
  unsigned int thread_count;
};

constexpr CudaEnumOpts default_opts = {50, .5, 4, 8, 32 * 256};

void search_enumeration_cuda(const double *mu, const double *rdiag,
                             const unsigned int enum_dimensions,
                             const float *start_point_coefficients, unsigned int start_point_count,
                             unsigned int start_point_dim, EvaluatorStrategy evaluator,
                             double initial_radius, CudaEnumOpts opts = default_opts);

#endif