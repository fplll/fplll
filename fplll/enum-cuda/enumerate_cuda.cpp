#include "enumerate_cuda.h"
#include <fplll/enum/enumerate_dyn.h>
#include "cuda_wrapper.h"
#include <iostream>

FPLLL_BEGIN_NAMESPACE

template<typename FT>
struct AllPointsEvaluator : public Evaluator<FT> {
    
  using Evaluator<FT>::max_sols;
  using Evaluator<FT>::strategy;
  using Evaluator<FT>::findsubsols;
  using Evaluator<FT>::normExp;
  using Evaluator<FT>::sol_count;
  using Evaluator<FT>::solutions;
  using Evaluator<FT>::size;

  AllPointsEvaluator()
      : Evaluator<FT>(std::numeric_limits<size_t>::max(), EvaluatorStrategy::EVALSTRATEGY_BEST_N_SOLUTIONS, false)
  {
  }
  virtual ~AllPointsEvaluator() {}

  virtual void eval_sol(const vector<FT> &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist) 
  {
    ++sol_count;
    solutions.emplace(new_partial_dist, new_sol_coord);
  }
  
  virtual void eval_sub_sol(int offset, const vector<FT> &new_sub_sol_coord,
                            const enumf &sub_dist) 
  {}

  typename Evaluator<FT>::container_t::const_iterator shortest_first_begin() const {
      return solutions.begin(); 
  }

  typename Evaluator<FT>::container_t::const_iterator shortest_first_end() const {
      return solutions.end();
  }
};

template <typename ZT, typename FT>
std::unique_ptr<AllPointsEvaluator<FT>> CudaEnumeration<ZT, FT>::enumerate_start_points(size_t start_dims, FT &fmaxdist, long fmaxdistexpo) {

    std::unique_ptr<AllPointsEvaluator<FT>> start_point_enumeration_evaluator(new AllPointsEvaluator<FT>());
    EnumerationDyn<ZT, FT> start_point_enumeration(_gso, *start_point_enumeration_evaluator);
    std::vector<FT> target_coord;
    std::vector<enumxt> subtree;
    std::vector<enumxt> pruning;
    start_point_enumeration.enumerate(_first + _d - start_dims, _first + _d, fmaxdist, fmaxdistexpo, target_coord, subtree, pruning, false, false); 
    return start_point_enumeration_evaluator;
}

template <typename ZT, typename FT>
bool CudaEnumeration<ZT, FT>::enumerate(int first, int last, FT &fmaxdist, long fmaxdistexpo,
                                        const vector<enumf> &pruning, bool dual)
{
    try {

        std::cout << "cuda enumeration called" << std::endl;

        if (dual) {
            throw "cuda enumeration does not yet support dual lattice enumeration";
        }
        _dual = false;

        if (!pruning.empty()) {
            throw "cuda enumeration does not yet support pruning";
        }
        _pruning = pruning;

        if (last == -1) {
            last = _gso.d;
        }
        _first = first;
        _d     = last - _first;
        
        cuenum::EvaluatorStrategy evaluator;
        FastErrorBoundedEvaluator* concrete_evaluator = dynamic_cast<FastErrorBoundedEvaluator*>(&_evaluator);
        if (concrete_evaluator != nullptr) {
            if (concrete_evaluator->max_sols > 1) {
                throw "cuda enumeration currently supports finding only one solution";
            } else if (concrete_evaluator->strategy != EvaluatorStrategy::EVALSTRATEGY_BEST_N_SOLUTIONS && concrete_evaluator->strategy != EvaluatorStrategy::EVALSTRATEGY_OPPORTUNISTIC_N_SOLUTIONS) {
                throw "cuda enumeration currently supports only finding the best solution";
            }
            evaluator = cuenum::EvaluatorStrategy::FAST;
        } else {
            throw "cuda enumeration does not accept arbitrary evaluators; using the cuda enumeration with a custom evaluator requires compiling enum.cuh with the correct template arguments using a cuda compiler";
        }

        std::unique_ptr<enumf[]> mu(new enumf[_d * _d]);
        std::unique_ptr<enumf[]> rdiag(new enumf[_d]);
        std::unique_ptr<enumf[]> pruning(new enumf[_d]);

        callback_set_config(mu.get(), _d, true, rdiag.get(), pruning.get());

        cuenum::CudaEnumOpts options = cuenum::default_opts;

        // currently, the enumerated dimension count must be divisible by dimensions_per_level, so
        // we have to adapt start_dims accordingly
        size_t start_dims = 5;
        while (start_dims % options.dimensions_per_level != 0) {
            ++start_dims;
        }
        if (_d < 0 || start_dims >= static_cast<size_t>(_d)) {
            // use fallback, as cuda enumeration in such small dimensions is too much overhead
            return false;
        }

        std::unique_ptr<AllPointsEvaluator<FT>> start_point_evaluator = enumerate_start_points(start_dims, fmaxdist, fmaxdistexpo);
        typedef typename AllPointsEvaluator<FT>::container_t::const_iterator start_point_iter;
        start_point_iter begin = start_point_evaluator->shortest_first_begin();
        start_point_iter end = start_point_evaluator->shortest_first_end();
        auto start_points = cuenum::create_start_point_array<start_point_iter>(start_point_evaluator->size(), start_dims, begin, end);

        cuenum::search_enumeration_cuda(mu.get(), rdiag.get(), _d - start_dims, start_points.get(), start_point_evaluator->size(), start_dims, evaluator, fmaxdist.get_d(), options);

        return true;

    } catch (const char*& error) {
        std::cerr << "Error in cuda enumeration, using fallback: " << error << std::endl;
    }
    return false;
}

template class CudaEnumeration<Z_NR<mpz_t>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class CudaEnumeration<Z_NR<mpz_t>, FP_NR<long double>>;
#endif

#ifdef FPLLL_WITH_QD
template class CudaEnumeration<Z_NR<mpz_t>, FP_NR<dd_real>>;

template class CudaEnumeration<Z_NR<mpz_t>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class CudaEnumeration<Z_NR<mpz_t>, FP_NR<dpe_t>>;
#endif

template class CudaEnumeration<Z_NR<mpz_t>, FP_NR<mpfr_t>>;

template class CudaEnumeration<Z_NR<long>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class CudaEnumeration<Z_NR<long>, FP_NR<long double>>;
#endif

#ifdef FPLLL_WITH_QD
template class CudaEnumeration<Z_NR<long>, FP_NR<dd_real>>;

template class CudaEnumeration<Z_NR<long>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class CudaEnumeration<Z_NR<long>, FP_NR<dpe_t>>;
#endif

template class CudaEnumeration<Z_NR<long>, FP_NR<mpfr_t>>;

FPLLL_END_NAMESPACE