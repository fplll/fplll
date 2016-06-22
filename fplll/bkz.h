/* Copyright (C) 2011 Xavier Pujol
   (C) 2014 Martin Albrecht.

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

#ifndef FPLLL_BKZ_H
#define FPLLL_BKZ_H

#include "bkz_params.h"
#include "enum/evaluator.h"
#include "enum/enumerate.h"
#include "lll.h"

FPLLL_BEGIN_NAMESPACE

template <class FT> class BKZAutoAbort
{
public:
  BKZAutoAbort(MatGSO<Integer, FT> &m, int num_rows, int start_row = 0)
      : m(m), old_slope(numeric_limits<double>::max()), no_dec(-1), num_rows(num_rows),
        start_row(start_row)
  {
  }
  bool test_abort(double scale = 1.0, int maxNoDec = 5);

private:
  MatGSO<Integer, FT> &m;
  double old_slope;
  int no_dec;
  int num_rows;
  int start_row;
};

/** Finds the slope of the curve fitted to the lengths of the vectors from
    startRow to stopRow. The slope gives an indication of the quality of the
    LLL-reduced basis.
*/
template <class FT> double get_current_slope(MatGSO<Integer, FT> &m, int startRow, int stopRow);

/**
   @brief Use the Gaussian Heuristic to compute a bound on the length
   of the shortest vector.

   @param max_dist         output
   @param max_dist_expo    exponent of output
   @param block_size       block size
   @param root_det
   @param gh_factor
   @return
*/

template <class FT>
void compute_gaussian_heuristic(FT &max_dist, long max_dist_expo, int block_size, const FT &root_det, double gh_factor);

/**
 * Compute the (squared) root determinant of the basis.
 */
template<class FT>
FT get_root_det(MatGSO<Integer, FT>& m, int start, int end);

/* The matrix must be LLL-reduced */
template <class FT> class BKZReduction
{
public:
  BKZReduction(MatGSO<Integer, FT> &m, LLLReduction<Integer, FT> &lllObj, const BKZParam &param);
  ~BKZReduction();

  bool svp_preprocessing(int kappa, int block_size, const BKZParam &param);

  bool svp_postprocessing(int kappa, int blockSize, const vector<FT> &solution);
  
  bool dsvp_postprocessing(int kappa, int block_size, const vector<FT> &solution);

  /**
     Run enumeration to find a new shortest vector in the sublattice B[kappa,kappa+blockSize]

     @param kappa Start row
     @param block_size Block size to use, this may be < param.blockSize
     @param param Parameters to use for this enumeration (blockSize is ignored)

     ``_ex`` variant is exception handling.
  */

  bool svp_reduction(int kappa, int block_size, const BKZParam &param);

  bool svp_reduction_Ex(int kappa, int block_size, const BKZParam &param, bool &clean)
  {
    try
    {
      clean = svp_reduction(kappa, block_size, param);
      return true;
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }
  }
  
  bool dsvp_reduction(int kappa, int block_size, const BKZParam &par);
  
  bool dsvp_reduction_Ex(int kappa, int block_size, const BKZParam &param, bool &clean)
  {
    try
    {
      clean = dsvp_reduction(kappa, block_size, param);
      return true;
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }
  }

  bool tour(const int loop, int &kappaMax, const BKZParam &param, int minRow, int maxRow);

  bool tour_ex(const int loop, int &kappaMax, const BKZParam &param, int minRow, int maxRow,
               bool &clean)
  {
    try
    {
      clean = tour(loop, kappaMax, param, minRow, maxRow);
      return true;
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }
  }
  
  bool sd_tour(const int loop, const BKZParam &param, int minRow, int maxRow);

  bool sd_tour_ex(const int loop, const BKZParam &param, int minRow, int maxRow,
               bool &clean)
  {
    try
    {
      clean = sd_tour(loop, param, minRow, maxRow);
      return true;
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }
  }
  
  bool hkz(int &kappaMax, const BKZParam &param, int min_row, int max_row);

  bool bkz();

  void rerandomize_block(int minRow, int maxRow, int density);

  /** I/O **/

  void dump_gso(const std::string filename, const std::string prefix, bool append = true);

  int status;

  /**
      Number of nodes visited during enumeration.
  */

  long nodes;

private:
  void print_tour(const int loop, int minRow, int maxRow);
  void print_params(const BKZParam &param, ostream &out);

  bool set_status(int newStatus);

  const Pruning &get_pruning(int kappa, int blockSize, const BKZParam &par) const;

  bool trunc_tour(int &kappaMax, const BKZParam &param, int minRow, int maxRow);
  bool trunc_dtour(const BKZParam &param, int minRow, int maxRow);

  const BKZParam &param;
  int num_rows;
  MatGSO<Integer, FT> &m;
  LLLReduction<Integer, FT> &lll_obj;
  FastEvaluator<FT> evaluator;
  FT delta;

  // Temporary data
  const vector<FT> empty_target, empty_sub_tree;
  FT max_dist, delta_max_dist;
  double cputime_start;
};

template <class FT>
bool BKZReduction<FT>::svp_preprocessing(int kappa, int block_size, const BKZParam &param)
{
  bool clean = true;

  FPLLL_DEBUG_CHECK(param.strategies.size() > block_size);

  auto &preproc = param.strategies[block_size].preprocessing_blocksizes;
  for (auto it = preproc.begin(); it != preproc.end(); ++it)
  {
    int dummy_kappa_max = num_rows;
    BKZParam prepar   = BKZParam(*it, param.strategies);
    clean &= tour(0, dummy_kappa_max, prepar, kappa, kappa + block_size);
  }

  return clean;
}

template <class FT>
bool BKZReduction<FT>::dsvp_postprocessing(int kappa, int block_size, const vector<FT> &solution)
{ 
  vector<FT> x = solution;
  
  int d = block_size;
  m.rowOpBegin(kappa, kappa + d);
  // don't want to deal with negativ coefficients
  for (int i = 0; i < d; i++) 
  {
    if (x[i] < 0) 
    {
      x[i].neg(x[i]);
      for (int j = 0; j < m.b.getCols(); j++) 
      {
        m.b[i + kappa][j].neg(m.b[i + kappa][j]);
      }
    }
  }
  
  // tree based gcd computation on x, performing operations also on b
  int off = 1;
  int k;
  while (off < d) 
  {
    k = d - 1;
    while(k - off >= 0) 
    {
      if (!(x[k].is_zero() && x[k - off].is_zero())) 
      {
        if (x[k] < x[k - off]) 
        {
          x[k].swap(x[k - off]);
          m.b.swapRows(kappa + k, kappa + k - off);
        }
        
        while (!x[k - off].is_zero()) 
        {
          while (x[k - off] <= x[k]) 
          {
            x[k] = x[k] - x[k - off];
            m.b[kappa + k].sub(m.b[kappa + k - off]);
          }
          
          x[k].swap(x[k - off]);
          m.b.swapRows(kappa + k, kappa + k - off);
        }
      }
      k -= 2 * off;
    }
    off *= 2;
  }
  
  m.rowOpEnd(kappa, kappa + d);
  if (!lll_obj.lll(kappa, kappa, kappa + d)) {
    return set_status(lll_obj.status);
  }
  return false;
}

template <class FT>
bool BKZReduction<FT>::dsvp_reduction(int kappa, int block_size, const BKZParam &par)
{
  bool clean = true;

  int lll_start = (par.flags & BKZ_BOUNDED_LLL) ? kappa : 0;

  if (!lll_obj.lll(lll_start, kappa, kappa + block_size))
  {
    throw lll_obj.status;
  }

  clean &= (lll_obj.nSwaps == 0);

  size_t trial                 = 0;
  double remaining_probability = 1.0;

  while (remaining_probability > 1. - par.min_success_probability)
  {
    if (trial > 0)
    {
      rerandomize_block(kappa, kappa + block_size - 1, par.rerandomization_density);
    }

    clean &= svp_preprocessing(kappa, block_size, par);

    long max_dist_expo;
    FT max_dist = m.getRExp(kappa + block_size - 1, kappa + block_size - 1, max_dist_expo);
    max_dist.pow_si(max_dist, -1);
    max_dist_expo *= -1;
    FT delta_max_dist;
    delta_max_dist = delta * max_dist;

    if ((par.flags & BKZ_GH_BND) && block_size > 30)
    {
      FT root_det = get_root_det(m, kappa, kappa + block_size);
      root_det.pow_si(root_det, -1);
      compute_gaussian_heuristic(max_dist, max_dist_expo, block_size, root_det, par.gh_factor);
    }

    const Pruning &pruning = get_pruning(kappa, block_size, par);

    vector<FT> &solCoord = evaluator.solCoord;
    solCoord.clear();
    Enumeration<FT> Enum(m, evaluator);
    Enum.enumerate( kappa, kappa + block_size, max_dist, max_dist_expo, vector<FT>(), vector<enumxt>(),
                    pruning.coefficients, true);
    nodes += Enum.getNodes();

    if (solCoord.empty())
    {
      if (pruning.coefficients[0] == 1 && !(par.flags & BKZ_GH_BND))
      {
        throw RED_ENUM_FAILURE;
      }
    }

    if (max_dist < delta_max_dist)
    {
      clean &= dsvp_postprocessing(kappa, block_size, solCoord);
    }
    remaining_probability *= (1 - pruning.probability);
    trial += 1;
  }
  return clean;
}

FPLLL_END_NAMESPACE

#endif
