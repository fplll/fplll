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

  bool svp_reduction(int kappa, int block_size, const BKZParam &param, bool dual = false);

  bool svp_reduction_Ex(int kappa, int block_size, const BKZParam &param, bool &clean, bool dual = false)
  {
    try
    {
      clean = svp_reduction(kappa, block_size, param, dual);
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

FPLLL_END_NAMESPACE

#endif
