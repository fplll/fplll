/* Copyright (C) 2008-2011 Xavier Pujol.

   This file is part of fplll. fplll is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   fplll is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with fplll. If not, see <http://www.gnu.org/licenses/>. */

#ifndef FPLLL_EVALUATOR_H
#define FPLLL_EVALUATOR_H

#include <cassert>
#include <fplll/gso_interface.h>
#include <fplll/util.h>
#include <functional>
#include <map>
#include <queue>

FPLLL_BEGIN_NAMESPACE

enum EvaluatorMode
{
  EVALMODE_SV    = 0,
  EVALMODE_CV    = 0,
  EVALMODE_COUNT = 1,
  EVALMODE_PRINT = 2
};

enum EvaluatorStrategy
{
  /*
   * Strategies for updating the enumeration bound as solutions are found.
   * Possible values are:
   *
   * EVALSTRATEGY_BEST_N_SOLUTIONS
   *   Starting with the max_sols-th solution, every time a new solution is found
   *   the enumeration bound is updated to the length of the longest solution. If
   *   more than max_sols were found, the longest is dropped.
   * EVALSTRATEGY_OPPORTUNISTIC_N_SOLUTIONS
   *   Every time a solution is found, update the enumeration distance to the length
   *   of the solution. If more than max_sols were found, the longest is dropped.
   * EVALSTRATEGY_FIRST_N_SOLUTIONS
   *   The enumeration bound is not updated. As soon as max_sols are found, enumeration
   *   stops.
   */
  EVALSTRATEGY_BEST_N_SOLUTIONS          = 0,
  EVALSTRATEGY_OPPORTUNISTIC_N_SOLUTIONS = 1,
  EVALSTRATEGY_FIRST_N_SOLUTIONS         = 2
};

/**
 * Evaluator stores the solutions found by enumerate, and updates the enumeration bound
 * It thus provides an interface to the enumerator,
 * as well as a basic interface to return solutions.
 * Specializations will implement specific behaviour and additional interfaces.
 */

template <class FT> class Evaluator
{
public:
  Evaluator(size_t nr_solutions               = 1,
            EvaluatorStrategy update_strategy = EVALSTRATEGY_BEST_N_SOLUTIONS,
            bool find_subsolutions            = false)
      : max_sols(nr_solutions), strategy(update_strategy), findsubsols(find_subsolutions),
        sol_count(0)
  {
    FPLLL_CHECK(nr_solutions > 0, "Evaluator: nr_solutions must be strictly positive!");
    FPLLL_CHECK(strategy <= 2, "Evaluator: invalid strategy");
  }
  virtual ~Evaluator() {}

  /** configuration */
  size_t max_sols;
  EvaluatorStrategy strategy;
  bool findsubsols;

  /** Solutions found in the lattice */
  // multimap storing solutions mapped by their length
  // the longest solution is the first, the shortest solution last
  typedef std::multimap<FT, std::vector<FT>, std::greater<FT>> container_t;
  container_t solutions;
  size_t sol_count;

  /** Subsolutions found in the lattice */
  std::vector<std::pair<FT, std::vector<FT>>> sub_solutions;

  /** interface to resulting solutions */
  typename container_t::const_reverse_iterator begin() const { return solutions.rbegin(); }
  typename container_t::const_reverse_iterator end() const { return solutions.rend(); }
  typename container_t::reverse_iterator begin() { return solutions.rbegin(); }
  typename container_t::reverse_iterator end() { return solutions.rend(); }
  size_t size() const { return solutions.size(); }
  bool empty() const { return solutions.empty(); }

  /** interface for the enumerator */
  virtual void eval_sol(const vector<FT> &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist) = 0;

  virtual void eval_sub_sol(int offset, const vector<FT> &new_sub_sol_coord,
                            const enumf &sub_dist) = 0;

  virtual void set_normexp(long norm_exp) { normExp = norm_exp; }
  long normExp;

protected:
  /** calculate enumeration bound based on dist */
  virtual enumf calc_enum_bound(const FT &dist) const
  {
    FT tmp;
    tmp.mul_2si(dist, -normExp);
    return tmp.get_d(GMP_RNDU);
  }

  /** processes solution into multimap and adjusts max_dist according to strategy */
  void process_sol(const FT &dist, const vector<FT> &coord, enumf &max_dist)
  {
    ++sol_count;
    solutions.emplace(dist, coord);
    switch (strategy)
    {
    case EVALSTRATEGY_BEST_N_SOLUTIONS:
      if (solutions.size() < max_sols)
        return;
      // remove the longest solution, and use the new longest dist to update max_dist
      if (solutions.size() > max_sols)
        solutions.erase(solutions.begin());
      max_dist = calc_enum_bound(solutions.begin()->first);
      break;

    case EVALSTRATEGY_OPPORTUNISTIC_N_SOLUTIONS:
      // always use dist to update max_dist
      max_dist = calc_enum_bound(dist);
      if (solutions.size() <= max_sols)
        return;
      // remove longest solution
      solutions.erase(solutions.begin());
      break;

    case EVALSTRATEGY_FIRST_N_SOLUTIONS:
      if (solutions.size() < max_sols)
        return;
      // when desired nr of solutions has been reached, set enum bound to zero
      max_dist = 0;
      break;

    default:
      FPLLL_CHECK(false, "Evaluator: invalid strategy switch!");
    }
  }
};

/**
 * Simple solution evaluator which provides a result without error bound.
 * The same instance can be used for several calls to enumerate on different
 * problems.
 */
template <class FT> class FastEvaluator : public Evaluator<FT>
{
public:
  using Evaluator<FT>::max_sols;
  using Evaluator<FT>::strategy;
  using Evaluator<FT>::findsubsols;
  using Evaluator<FT>::normExp;
  using Evaluator<FT>::sub_solutions;

  FastEvaluator(size_t nr_solutions               = 1,
                EvaluatorStrategy update_strategy = EVALSTRATEGY_BEST_N_SOLUTIONS,
                bool find_subsolutions            = false)
      : Evaluator<FT>(nr_solutions, update_strategy, find_subsolutions)
  {
  }
  virtual ~FastEvaluator() {}

  virtual void eval_sol(const vector<FT> &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist)
  {
    FT dist = new_partial_dist;
    dist.mul_2si(dist, normExp);

    // store solution and update max_dist according to strategy
    this->process_sol(dist, new_sol_coord, max_dist);
  }

  virtual void eval_sub_sol(int offset, const vector<FT> &new_sub_sol_coord, const enumf &sub_dist)
  {
    FT dist = sub_dist;
    dist.mul_2si(dist, normExp);

    sub_solutions.resize(std::max(sub_solutions.size(), std::size_t(offset + 1)));

    if (sub_solutions[offset].second.empty() || dist < sub_solutions[offset].first)
    {
      sub_solutions[offset].first  = dist;
      sub_solutions[offset].second = new_sub_sol_coord;
      for (int i = 0; i < offset; ++i)
        sub_solutions[offset].second[i] = 0.0;
    }
  }
};

/**
   @brief Callback function used by CallbackEvaluator.

 */

typedef bool(callback_evaluator_callback)(size_t n, enumf *new_sol_coord, void *ctx);

/**
   @brief A FastEvaluator which additionally checks whether the predicate ``callbackf(solution,
   ctx)`` accepts or rejects.

   @example tests/test_enum.cpp

 */
template <class FT> class CallbackEvaluator : public FastEvaluator<FT>
{

  std::function<callback_evaluator_callback> callbackf;
  void *ctx;
  enumf new_sol_coordf[FPLLL_MAX_ENUM_DIMENSION];

public:
  using FastEvaluator<FT>::max_sols;
  using FastEvaluator<FT>::strategy;
  using FastEvaluator<FT>::findsubsols;
  using FastEvaluator<FT>::normExp;
  using FastEvaluator<FT>::sub_solutions;

  CallbackEvaluator(std::function<callback_evaluator_callback> callbackf, void *ctx = NULL,
                    size_t nr_solutions               = 1,
                    EvaluatorStrategy update_strategy = EVALSTRATEGY_BEST_N_SOLUTIONS,
                    bool find_subsolutions            = false

                    )
      : FastEvaluator<FT>(nr_solutions, update_strategy, find_subsolutions), callbackf(callbackf),
        ctx(ctx)
  {
  }
  virtual ~CallbackEvaluator() {}

  virtual void eval_sol(const vector<FT> &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist)
  {
    assert(new_sol_coord.size() <= FPLLL_MAX_ENUM_DIMENSION);
    for (size_t i = 0; i < new_sol_coord.size(); i++)
    {
      new_sol_coordf[i] = new_sol_coord[i].get_d();
    }
    if (!callbackf(new_sol_coord.size(), new_sol_coordf, this->ctx))
      return;

    FastEvaluator<FT>::eval_sol(new_sol_coord, new_partial_dist, max_dist);
  }
};

/**
 * ErrorBoundEvaluator provides an extra interface to provide
 * information about the accuracy of solutions.
 */
class ErrorBoundedEvaluator : public Evaluator<FP_NR<mpfr_t>>
{
public:
  ErrorBoundedEvaluator(int dim, const Matrix<FP_NR<mpfr_t>> &mmu, const Matrix<FP_NR<mpfr_t>> &mr,
                        EvaluatorMode evalmode, size_t nr_solutions = 1,
                        EvaluatorStrategy update_strategy = EVALSTRATEGY_BEST_N_SOLUTIONS,
                        bool find_subsolutions            = false)
      : Evaluator(nr_solutions, update_strategy, find_subsolutions), eval_mode(evalmode), d(dim),
        mu(mmu), r(mr), input_error_defined(false)
  {
    max_dr_diag.resize(d);
    max_dm_u.resize(d);
  }

  virtual ~ErrorBoundedEvaluator() {}

  /** Configuration */
  EvaluatorMode eval_mode;
  int d;
  const Matrix<FP_NR<mpfr_t>> &mu;
  const Matrix<FP_NR<mpfr_t>> &r;

  /* To enable error estimation, the caller must set
  input_error_defined=true and fill max_dr_diag and max_dm_u */
  bool input_error_defined;
  vector<FP_NR<mpfr_t>> max_dr_diag, max_dm_u;  // Error bounds on input parameters
  //  FP_NR<mpfr_t> last_partial_dist;          // Approx. squared norm of the last solution

  void init_delta_def(int prec, double rho, bool withRoundingToEnumf);

  /**
   * Computes max_error such that
   * normOfSolution^2 <= (1 + max_error) * lambda_1(L)^2.
   * The default implementation might fail (i.e. return false).
   */
  virtual bool get_max_error(FP_NR<mpfr_t> &max_error, const FP_NR<mpfr_t> &sol_dist) = 0;

  // Internal use
  bool get_max_error_aux(const FP_NR<mpfr_t> &max_dist, bool boundOnExactVal, FP_NR<mpfr_t> &maxDE);
};

/**
 * Simple solution evaluator which provides a non-certified result, but can
 * give an error bound.
 * The same object can be used for several calls to enumerate on different
 * instances.
 */
class FastErrorBoundedEvaluator : public ErrorBoundedEvaluator
{
public:
  FastErrorBoundedEvaluator(int d = 0, const Matrix<FP_NR<mpfr_t>> &mu = Matrix<FP_NR<mpfr_t>>(),
                            const Matrix<FP_NR<mpfr_t>> &r = Matrix<FP_NR<mpfr_t>>(),
                            EvaluatorMode eval_mode = EVALMODE_SV, size_t nr_solutions = 1,
                            EvaluatorStrategy update_strategy = EVALSTRATEGY_BEST_N_SOLUTIONS,
                            bool find_subsolutions            = false)
      : ErrorBoundedEvaluator(d, mu, r, eval_mode, nr_solutions, update_strategy, find_subsolutions)
  {
  }
  virtual ~FastErrorBoundedEvaluator() {}

  virtual bool get_max_error(FP_NR<mpfr_t> &max_error, const FP_NR<mpfr_t> &sol_dist);
  virtual void eval_sol(const vector<FP_NR<mpfr_t>> &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist);
  virtual void eval_sub_sol(int offset, const vector<FP_NR<mpfr_t>> &new_sub_sol_coord,
                            const enumf &sub_dist);
};

/**
 * ExactEvaluator stores the best solution found by enumerate.
 * The result is guaranteed, but the the evaluation of new solutions is longer.
 */
class ExactErrorBoundedEvaluator : public ErrorBoundedEvaluator
{
public:
  ExactErrorBoundedEvaluator(MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &_gso,
                             EvaluatorMode eval_mode, size_t nr_solutions = 1,
                             EvaluatorStrategy update_strategy = EVALSTRATEGY_BEST_N_SOLUTIONS,
                             bool find_subsolutions            = false)
      : ErrorBoundedEvaluator(_gso.d, _gso.get_mu_matrix(), _gso.get_r_matrix(), eval_mode,
                              nr_solutions, update_strategy, find_subsolutions),
        gso(_gso)
  {
    int_max_dist = -1;
  }

  virtual ~ExactErrorBoundedEvaluator() {}

  /**
   * Sets max_error to 0: the result is guaranteed.
   */
  virtual bool get_max_error(FP_NR<mpfr_t> &max_error, const FP_NR<mpfr_t> &sol_dist);

  virtual void eval_sol(const vector<FP_NR<mpfr_t>> &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist);

  virtual void eval_sub_sol(int offset, const vector<FP_NR<mpfr_t>> &new_sub_sol_coord,
                            const enumf &sub_dist);

  Z_NR<mpz_t> int_max_dist;  // Exact norm of the last vector

  Z_NR<mpz_t> exact_sol_dist(const vector<FP_NR<mpfr_t>> &sol_coord);
  Z_NR<mpz_t> exact_subsol_dist(int offset, const vector<FP_NR<mpfr_t>> &sol_coord);

private:
  FP_NR<mpfr_t> int_dist2Float(Z_NR<mpz_t> int_dist);
  MatGSOInterface<Z_NR<mpz_t>, FP_NR<mpfr_t>> &gso;
  // const ZZ_mat<mpz_t> &matrix;  // matrix of the lattice
};

FPLLL_END_NAMESPACE

#endif
