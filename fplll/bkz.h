/* Copyright (C) 2011 Xavier Pujol
   (C) 2014-2016 Martin R. Albrecht
   (C) 2016 Michael Walter

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

#include "bkz_param.h"
#include "enum/enumerate.h"
#include "enum/evaluator.h"
#include "lll.h"

FPLLL_BEGIN_NAMESPACE

template <class FT> class BKZAutoAbort
{
public:
  /**
     @brief

     @param m
     @param num_rows
     @param start_row
     @return
  */

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

/* The matrix must be LLL-reduced */
template <class FT> class BKZReduction
{
public:
  BKZReduction(MatGSO<Integer, FT> &m, LLLReduction<Integer, FT> &lll_obj, const BKZParam &param);
  ~BKZReduction();

  bool svp_preprocessing(int kappa, int block_size, const BKZParam &param);

  bool svp_postprocessing(int kappa, int block_size, const vector<FT> &solution, bool dual = false);
  
  bool svp_post_general(int kappa, int block_size, const vector<FT> &solution, bool dual = false);

  /**
     Run enumeration to find a new shortest vector in the sublattice B[kappa,kappa+block_size]

     @param kappa Start row
     @param block_size Block size to use, this may be < param.block_size
     @param param Parameters to use for this enumeration (block_size is ignored)

     ``_ex`` variant is exception handling.
  */

  bool svp_reduction(int kappa, int block_size, const BKZParam &param, bool dual = false);

  bool svp_reduction_ex(int kappa, int block_size, const BKZParam &param, bool &clean,
                        bool dual = false)
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

  bool tour(const int loop, int &kappa_max, const BKZParam &param, int min_row, int max_row);

  bool tour_ex(const int loop, int &kappa_max, const BKZParam &param, int min_row, int max_row,
               bool &clean)
  {
    try
    {
      clean = tour(loop, kappa_max, param, min_row, max_row);
      return true;
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }
  }

  bool sd_tour(const int loop, const BKZParam &param, int min_row, int max_row);

  bool sd_tour_ex(const int loop, const BKZParam &param, int min_row, int max_row, bool &clean)
  {
    try
    {
      clean = sd_tour(loop, param, min_row, max_row);
      return true;
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }
  }

  bool hkz(int &kappa_max, const BKZParam &param, int min_row, int max_row);

  bool hkz_ex(int &kappa_max, const BKZParam &param, int min_row, int max_row, bool &clean)
  {
    try
    {
      clean = hkz(kappa_max, param, min_row, max_row);
      return true;
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }
  }

  bool slide_tour(const int loop, const BKZParam &param, int min_row, int max_row);

  bool slide_tour_ex(const int loop, const BKZParam &param, int min_row, int max_row, bool &clean)
  {
    try
    {
      clean = slide_tour(loop, param, min_row, max_row);
      return true;
    }
    catch (RedStatus &e)
    {
      return set_status(e);
    }
  }

  bool bkz();

  /** Randomize basis between from ``min_row`` and ``max_row`` (exclusive)

      1. permute rows
      2. apply lower triangular matrix with coefficients in -1,0,1
      3. LLL reduce result

      @param min_row start in this row

      @param max_row stop at this row (exclusive)

      @param density number of non-zero coefficients in lower triangular
      transformation matrix
  **/

  void rerandomize_block(int min_row, int max_row, int density);

  /** I/O **/

  void dump_gso(const std::string &filename, const std::string &prefix, bool append = true);

  int status;

  /**
      Number of nodes visited during enumeration.
  */

  long nodes;

private:
  void print_tour(const int loop, int min_row, int max_row);
  void print_params(const BKZParam &param, ostream &out);

  bool set_status(int new_status);

  const Pruning &get_pruning(int kappa, int block_size, const BKZParam &par) const;

  bool trunc_tour(int &kappa_max, const BKZParam &param, int min_row, int max_row);
  bool trunc_dtour(const BKZParam &param, int min_row, int max_row);

  const BKZParam &param;
  int num_rows;
  MatGSO<Integer, FT> &m;
  LLLReduction<Integer, FT> &lll_obj;
  FastEvaluator<FT> evaluator;
  FT delta;

  const char *algorithm;
  // Temporary data
  const vector<FT> empty_target, empty_sub_tree;
  FT max_dist, delta_max_dist;
  double cputime_start;
  FT sld_potential;
  FT ftmp;
};

int bkz_reduction(IntMatrix *B, IntMatrix *U, const BKZParam &param,
                  FloatType float_type = FT_DEFAULT, int precision = 0);
int bkz_reduction(IntMatrix &b, int block_size, int flags = BKZ_DEFAULT,
                  FloatType float_type = FT_DEFAULT, int precision = 0);
int bkz_reduction(IntMatrix &b, IntMatrix &u, int block_size, int flags = BKZ_DEFAULT,
                  FloatType float_type = FT_DEFAULT, int precision = 0);

int hkz_reduction(IntMatrix &b, int flags = HKZ_DEFAULT, FloatType float_type = FT_DEFAULT,
                  int precision = 0);

FPLLL_END_NAMESPACE

#endif
