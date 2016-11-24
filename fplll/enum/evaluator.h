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

#include "../util.h"
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

/**
 * Evaluator stores the best solution found by enumerate. The Float
 * specialization provides additional information about the solution accuracy.
 */
template <class FT> class Evaluator
{
public:
  Evaluator(size_t max_aux_solutions = 0, bool find_subsolutions = false, bool always_update_radius = true)
      : max_aux_sols(max_aux_solutions), findsubsols(find_subsolutions), new_sol_flag(false)
  {
    always_update_rad = always_update_radius||max_aux_solutions==0;
  }
  virtual ~Evaluator() {}

  /** Called by enumerate when a solution is found.
     Input: new_sol_coord = coordinates of the solution in Gram-Schmidt basis
     new_partial_dist = estimated distance between the solution and the
     orthogonal projection of target on the lattice space
     max_dist = current bound of the algorithm
     Output: max_dist can be decreased */
  virtual void eval_sol(const vector<FT> &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist) = 0;

  virtual void eval_sub_sol(int offset, const vector<FT> &new_sub_sol_coord,
                            const enumf &sub_dist) = 0;

  virtual void set_normexp(long /*norm_exp*/) {}

  /** Coordinates of the solution in the lattice */
  vector<FT> sol_coord;
  enumf sol_dist;

  /** Other solutions found in the lattice */
  size_t max_aux_sols;
  bool always_update_rad;
  std::multimap<enumf, vector<FT>, std::greater<enumf> > aux_sols;
  virtual std::vector<std::pair<enumf, vector<FT> > > multimap2pairs()
  {
      std::vector<std::pair<enumf, vector<FT> > > result;
      for (auto it = aux_sols.rbegin(), itend = aux_sols.rend(); it != itend; ++it)
          result.emplace_back(it->first, it->second);
      return result;
  }
  virtual void set_max_aux_sols(size_t new_max)
  {
      FPLLL_CHECK(aux_sols.size() == 0, "invalid call");
      max_aux_sols = new_max;
      always_update_rad = always_update_rad||max_aux_sols==0;
  }

  /** Subsolutions found in the lattice */
  bool findsubsols;
  vector<vector<FT>> sub_sol_coord;
  vector<enumf> sub_sol_dist;

  /** Set to true when sol_coord is updated */
  bool new_sol_flag;
};

/**
 * Simple solution evaluator which provides a result without error bound.
 * The same instance can be used for several calls to enumerate on different
 * problems.
 */
template <class FT> class FastEvaluator : public Evaluator<FT>
{
public:
  using Evaluator<FT>::sol_coord;
  using Evaluator<FT>::sol_dist;
  using Evaluator<FT>::new_sol_flag;
  using Evaluator<FT>::aux_sols;
  using Evaluator<FT>::sub_sol_coord;
  using Evaluator<FT>::sub_sol_dist;
  using Evaluator<FT>::max_aux_sols;
  using Evaluator<FT>::always_update_rad;

  FastEvaluator(size_t max_aux_solutions = 0, bool find_subsolutions = false, bool always_update_radius = true)
      : Evaluator<FT>(max_aux_solutions, find_subsolutions, always_update_radius)
  {
  }

  virtual ~FastEvaluator() {}

  /**
   * Called by enumerate when a solution is found.
   * FastEvaluator always accepts the solution and sets the bound max_dist to
   * new_partial_dist.
   *
   * @param new_sol_coord    Coordinates of the solution in the lattice
   * @param new_partial_dist Floating-point evaluation of the norm of the solution
   * @param max_dist        Bound of the enumeration (updated by the function)
   * @param normExp        r(i, i) is divided by 2^normExp in enumerate before
   *                       being converted to double
   */
  virtual void eval_sol(const vector<FT> &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist)
  {
    if (max_aux_sols != 0 && !sol_coord.empty())
    {
      aux_sols.emplace(sol_dist, sol_coord);
      if (aux_sols.size() > max_aux_sols)
      {
        max_dist = aux_sols.erase(aux_sols.begin())->first;
      }
    }
    sol_coord = new_sol_coord;
    sol_dist = new_partial_dist;
    new_sol_flag        = true;
    if (always_update_rad)
    {
        max_dist = sol_dist;
    }
  }

  virtual void eval_sub_sol(int offset, const vector<FT> &new_sub_sol_coord, const enumf &sub_dist)
  {
    sub_sol_coord.resize(std::max(sub_sol_coord.size(), std::size_t(offset + 1)));
    sub_sol_dist.resize(sub_sol_coord.size(), -1.0);
    if (sub_sol_dist[offset] == -1.0 || sub_dist < sub_sol_dist[offset])
    {
      sub_sol_coord[offset] = new_sub_sol_coord;
      for (int i                 = 0; i < offset; ++i)
        sub_sol_coord[offset][i] = 0.0;
      sub_sol_dist[offset]       = sub_dist;
    }
  }
};

/**
 * Evaluator stores the best solution found by enumerate and provides
 * information about the accuracy of this solution.
 */
template <> class Evaluator<Float>
{
public:
  Evaluator<Float>(int d, const Matrix<Float> &mu, const Matrix<Float> &r, int eval_mode,
                   size_t max_aux_solutions = 0, bool find_subsolutions = false, bool always_update_radius = true)
      : max_aux_sols(max_aux_solutions), findsubsols(find_subsolutions), new_sol_flag(false),
        eval_mode(eval_mode), input_error_defined(false), d(d), mu(mu), r(r)
  {
    always_update_rad = always_update_radius||max_aux_solutions==0;
    max_dr_diag.resize(d);
    max_dm_u.resize(d);
  }

  virtual ~Evaluator<Float>() {}

  virtual void set_normexp(long norm_exp) { normExp = norm_exp; }
  long normExp;

  void init_delta_def(int prec, double rho, bool withRoundingToEnumf);

  /**
   * Computes max_error such that
   * normOfSolution^2 <= (1 + max_error) * lambda_1(L)^2.
   * The default implementation might fail (i.e. return false).
   */
  virtual bool get_max_error(Float &max_error) = 0;

  /**
   * Called by enumerate when a solution is found.
   * The default implementation always accepts the solution and sets the bound
   * max_dist to new_partial_dist.
   *
   * @param new_sol_coord    Coordinates of the solution
   * @param new_partial_dist Floating-point estimation of the norm of the solution
   * @param max_dist        Bound of the enumeration (updated by the function)
   * @param normExp        It is assumed that r(i, i) is divided by 2^normExp
   *                       in enumerate
   */
  virtual void eval_sol(const FloatVect &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist) = 0;
  virtual void eval_sub_sol(int offset, const FloatVect &new_sub_sol_coord,
                            const enumf &sub_dist) = 0;

  virtual std::vector<std::pair<enumf, vector<Float> > > multimap2pairs()
  {
      std::vector<std::pair<enumf, vector<Float> > > result;
      for (auto it = aux_sols.rbegin(), itend = aux_sols.rend(); it != itend; ++it)
          result.emplace_back(it->first, it->second);
      return result;
  }
  virtual void set_max_aux_sols(size_t new_max)
  {
      FPLLL_CHECK(aux_sols.size() == 0, "invalid call");
      max_aux_sols = new_max;
      always_update_rad = always_update_rad||max_aux_sols==0;
  }


  // Internal use
  bool get_max_error_aux(const Float &max_dist, bool boundOnExactVal, Float &maxDE);

  /** Coordinates of the solution in the lattice */
  FloatVect sol_coord;
  enumf sol_dist;

  /** Other solutions found in the lattice */
  size_t max_aux_sols;
  bool always_update_rad;
  std::multimap<enumf, FloatVect, std::greater<enumf> > aux_sols;

  /** Subsolutions found in the lattice */
  bool findsubsols;
  vector<FloatVect> sub_sol_coord;
  vector<enumf> sub_sol_dist;

  /** Set to true when sol_coord is updated */
  bool new_sol_flag;
  /** Incremented when sol_coord is updated */
  long long sol_count;
  int eval_mode;

  /* To enable error estimation, the caller must set
     input_error_defined=true and fill max_dr_diag and max_dm_u */
  bool input_error_defined;
  FloatVect max_dr_diag, max_dm_u;  // Error bounds on input parameters
  Float last_partial_dist;          // Approx. squared norm of the last solution

  int d;
  const Matrix<Float> &mu;
  const Matrix<Float> &r;
};

/**
 * Simple solution evaluator which provides a non-certified result, but can
 * give an error bound.
 * The same object can be used for several calls to enumerate on different
 * instances.
 */
template <> class FastEvaluator<Float> : public Evaluator<Float>
{
public:
  FastEvaluator(int d = 0, const Matrix<Float> &mu = Matrix<Float>(),
                const Matrix<Float> &r = Matrix<Float>(), int eval_mode = EVALMODE_SV,
                size_t max_aux_solutions = 0, bool find_subsolutions = false)
      : Evaluator<Float>(d, mu, r, eval_mode, max_aux_solutions, find_subsolutions)
  {
  }
  virtual ~FastEvaluator() {}

  virtual bool get_max_error(Float &max_error);
  virtual void eval_sol(const FloatVect &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist);
  virtual void eval_sub_sol(int offset, const FloatVect &new_sub_sol_coord, const enumf &sub_dist);
};

/**
 * ExactEvaluator stores the best solution found by enumerate.
 * The result is guaranteed, but the the evaluation of new solutions is longer.
 */
class ExactEvaluator : public Evaluator<Float>
{
public:
  ExactEvaluator(int d, const IntMatrix &matrix, const Matrix<Float> &mu, const Matrix<Float> &r,
                 int eval_mode, size_t max_aux_solutions = 0, bool find_subsolutions = false)
      : Evaluator<Float>(d, mu, r, eval_mode, max_aux_solutions, find_subsolutions), matrix(matrix)
  {
    int_max_dist = -1;
  }

  /**
   * Sets max_error to 0: the result is guaranteed.
   */
  virtual bool get_max_error(Float &max_error);

  virtual void eval_sol(const FloatVect &new_sol_coord, const enumf &new_partial_dist,
                        enumf &max_dist);

  virtual void eval_sub_sol(int offset, const FloatVect &new_sub_sol_coord, const enumf &sub_dist);

  Integer int_max_dist;  // Exact norm of the last vector

  std::priority_queue<Integer> aux_sol_int_dist;  // Exact norm of aux vectors
  vector<Integer> sub_sol_int_dist;      // Exact norm of sub vectors

private:
  enumf int_dist2enumf(Integer int_dist);

  const IntMatrix &matrix;  // matrix of the lattice
};

FPLLL_END_NAMESPACE

#endif
