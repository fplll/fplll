/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2011 Xavier Pujol.

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

/* Template source file */

#include "gso_interface.h"

FPLLL_BEGIN_NAMESPACE

template <class ZT, class FT>
inline void MatGSOInterface<ZT, FT>::invalidate_gso_row(int i, int new_valid_cols)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && new_valid_cols >= 0 && new_valid_cols <= i + 1);
  gso_valid_cols[i] = min(gso_valid_cols[i], new_valid_cols);
}

template <class ZT, class FT> void MatGSOInterface<ZT, FT>::row_op_end(int first, int last)
{
#ifdef DEBUG
  FPLLL_DEBUG_CHECK(row_op_first == first && row_op_last == last);
  row_op_first = row_op_last = -1;
#endif
  for (int i = first; i < last; i++)
  {
    if (!enable_int_gram)
    {
      update_bf(i);
      invalidate_gram_row(i);
      for (int j = i + 1; j < n_known_rows; j++)
        gf(j, i).set_nan();
    }
    invalidate_gso_row(i, 0);
  }
  for (int i = last; i < n_known_rows; i++)
  {
    invalidate_gso_row(i, first);
  }
}

template <class ZT, class FT> inline ZT MatGSOInterface<ZT, FT>::get_max_gram()
{
  ZT tmp;
  if (enable_int_gram)
  {
    if (gptr == nullptr)
    {
      throw std::runtime_error("Error: gptr is equal to the nullpointer.");
    }
    Matrix<ZT> gr = *gptr;
    tmp           = gr(0, 0);
    for (int i = 0; i < n_known_rows; i++)
      tmp = tmp.max_z(gr(i, i));
  }
  else
  {
    FT tmp1 = gf(0, 0);
    for (int i = 0; i < n_known_rows; i++)
      tmp1 = tmp1.max_f(gf(i, i));
    tmp.set_f(tmp1);
  }
  return tmp;
}

template <class ZT, class FT> inline FT MatGSOInterface<ZT, FT>::get_max_bstar()
{
  FT tmp;
  tmp = r(0, 0);
  for (int i = 0; i < n_known_rows; i++)
    tmp = tmp.max_f(r(i, i));
  return tmp;
}

template <class ZT, class FT> long MatGSOInterface<ZT, FT>::get_max_mu_exp(int i, int n_columns)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && gso_valid_cols[i] >= n_columns);
  long max_expo = LONG_MIN, expo;
  for (int j = 0; j < n_columns; j++)
  {
    long expo2 = get_mu_exp(i, j, expo).exponent();
    max_expo   = max(max_expo, expo + expo2);
  }
  return max_expo;
}
// This function *can* be included in the interface; it has no base-class specific functions.
// If this function is going to be included, it must be removed from gso.* and gsogram.*
// and also added to gso_interface.h as non-virtual.
/*
template <class ZT, class FT>
void MatGSOInterface<ZT, FT>::row_addmul_we(int i, int j, const FT &x, long expo_add)
{
  FPLLL_DEBUG_CHECK(j >= 0 && i < n_known_rows && j < n_source_rows);
  long expo;
  long lx = x.get_si_exp_we(expo, expo_add);

  if (expo == 0)
  {
    if (lx == 1)
      row_add(i, j);
    else if (lx == -1)
      row_sub(i, j);
    else if (lx != 0)
      row_addmul_si(i, j, lx);
  }
  else if (row_op_force_long)
  {
    row_addmul_si_2exp(i, j, lx, expo);
  }
  else
  {
    x.get_z_exp_we(ztmp2, expo, expo_add);
    row_addmul_2exp(i, j, ztmp2, expo);
  }
}
*/

template <class ZT, class FT> bool MatGSOInterface<ZT, FT>::update_gso_row(int i, int last_j)
{
  // FPLLL_TRACE_IN("Updating GSO up to (" << i << ", " << last_j << ")");
  // FPLLL_TRACE("n_known_rows=" << n_known_rows << " n_source_rows=" << n_source_rows);
  if (i >= n_known_rows)
  {
    discover_row();
  }
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && last_j >= 0 && last_j < n_source_rows);

  int j = max(0, gso_valid_cols[i]);

  for (; j <= last_j; j++)
  {
    get_gram(ftmp1, i, j);
    FPLLL_DEBUG_CHECK(j == i || gso_valid_cols[j] >= j);
    for (int k = 0; k < j; k++)
    {
      ftmp2.mul(mu(j, k), r(i, k));
      ftmp1.sub(ftmp1, ftmp2);
    }
    r(i, j) = ftmp1;
    if (i > j)
    {
      mu(i, j).div(ftmp1, r(j, j));
      if (!mu(i, j).is_finite())
        return false;
    }
  }

  gso_valid_cols[i] = j;  // = max(0, gso_valid_cols[i], last_j + 1)
  // FPLLL_TRACE_OUT("End of GSO update");
  return true;
}

template <class ZT, class FT> void MatGSOInterface<ZT, FT>::lock_cols() { cols_locked = true; }

template <class ZT, class FT> void MatGSOInterface<ZT, FT>::unlock_cols()
{
  n_known_rows = n_source_rows;
  cols_locked  = false;
}

template <class ZT, class FT>
void MatGSOInterface<ZT, FT>::apply_transform(const Matrix<FT> &transform, int src_base,
                                              int target_base)
{
  int target_size = transform.get_rows(), src_size = transform.get_cols();
  int old_d = d;
  create_rows(target_size);
  for (int i = 0; i < target_size; i++)
  {
    for (int j = 0; j < src_size; j++)
    {
      row_addmul(old_d + i, src_base + j, transform(i, j));
    }
  }
  row_op_begin(target_base, target_base + target_size);
  for (int i = 0; i < target_size; i++)
  {
    row_swap(target_base + i, old_d + i);
  }
  row_op_end(target_base, target_base + target_size);
  remove_last_rows(target_size);
}

template <class ZT, class FT>
double MatGSOInterface<ZT, FT>::get_current_slope(int start_row, int stop_row)
{
  FT f, log_f;
  long expo;
  vector<double> x;
  x.resize(stop_row);
  for (int i = start_row; i < stop_row; i++)
  {
    update_gso_row(i);
    f = get_r_exp(i, i, expo);
    log_f.log(f, GMP_RNDU);
    x[i] = log_f.get_d() + expo * std::log(2.0);
  }
  int n         = stop_row - start_row;
  double i_mean = (n - 1) * 0.5 + start_row, x_mean = 0, v1 = 0, v2 = 0;
  for (int i = start_row; i < stop_row; i++)
  {
    x_mean += x[i];
  }
  x_mean /= n;
  for (int i = start_row; i < stop_row; i++)
  {
    v1 += (i - i_mean) * (x[i] - x_mean);
    v2 += (i - i_mean) * (i - i_mean);
  }
  return v1 / v2;
}

template <class ZT, class FT> FT MatGSOInterface<ZT, FT>::get_root_det(int start_row, int end_row)
{
  start_row   = max(0, start_row);
  end_row     = min(d, end_row);
  FT h        = (double)(end_row - start_row);
  FT root_det = get_log_det(start_row, end_row) / h;
  root_det.exponential(root_det);
  return root_det;
}

template <class ZT, class FT> FT MatGSOInterface<ZT, FT>::get_log_det(int start_row, int end_row)
{
  FT log_det = 0.0;
  start_row  = max(0, start_row);
  end_row    = min(d, end_row);
  FT h;
  for (int i = start_row; i < end_row; ++i)
  {
    get_r(h, i, i);
    log_det += log(h);
  }
  return log_det;
}

template <class ZT, class FT>
FT MatGSOInterface<ZT, FT>::get_slide_potential(int start_row, int end_row, int block_size)
{
  FT potential = 0.0;
  int p        = (end_row - start_row) / block_size;
  if ((end_row - start_row) % block_size == 0)
  {
    --p;
  }
  for (int i = 0; i < p; ++i)
  {
    potential += (p - i) * get_log_det(i * block_size, (i + 1) * block_size);
  }
  return potential;
}

template <class FT>
void adjust_radius_to_gh_bound(FT &max_dist, long max_dist_expo, int block_size, const FT &root_det,
                               double gh_factor)
{
  double t = (double)block_size / 2.0 + 1;
  t        = lgamma(t);
  t        = pow(M_E, t * 2.0 / (double)block_size);
  t        = t / M_PI;
  FT f     = t;
  f        = f * root_det;
  f.mul_2si(f, -max_dist_expo);
  f = f * gh_factor;
  if (f < max_dist)
  {
    max_dist = f;
  }
}

template class MatGSOInterface<Z_NR<long>, FP_NR<double>>;
template class MatGSOInterface<Z_NR<double>, FP_NR<double>>;
template class MatGSOInterface<Z_NR<>, FP_NR<double>>;
template void adjust_radius_to_gh_bound<FP_NR<double>>(FP_NR<double> &max_dist, long max_dist_expo,
                                                       int block_size,
                                                       const FP_NR<double> &root_det,
                                                       double gh_factor);

#ifdef FPLLL_WITH_LONG_DOUBLE
template class MatGSOInterface<Z_NR<long>, FP_NR<long double>>;
template class MatGSOInterface<Z_NR<double>, FP_NR<long double>>;
template class MatGSOInterface<Z_NR<>, FP_NR<long double>>;
template void adjust_radius_to_gh_bound<FP_NR<long double>>(FP_NR<long double> &max_dist,
                                                            long max_dist_expo, int block_size,
                                                            const FP_NR<long double> &root_det,
                                                            double gh_factor);

#endif

#ifdef FPLLL_WITH_QD
template class MatGSOInterface<Z_NR<long>, FP_NR<dd_real>>;
template class MatGSOInterface<Z_NR<double>, FP_NR<dd_real>>;
template class MatGSOInterface<Z_NR<>, FP_NR<dd_real>>;
template void adjust_radius_to_gh_bound<FP_NR<dd_real>>(FP_NR<dd_real> &max_dist,
                                                        long max_dist_expo, int block_size,
                                                        const FP_NR<dd_real> &root_det,
                                                        double gh_factor);
template class MatGSOInterface<Z_NR<long>, FP_NR<qd_real>>;
template class MatGSOInterface<Z_NR<double>, FP_NR<qd_real>>;
template class MatGSOInterface<Z_NR<>, FP_NR<qd_real>>;
template void adjust_radius_to_gh_bound<FP_NR<qd_real>>(FP_NR<qd_real> &max_dist,
                                                        long max_dist_expo, int block_size,
                                                        const FP_NR<qd_real> &root_det,
                                                        double gh_factor);
#endif

#ifdef FPLLL_WITH_DPE
template class MatGSOInterface<Z_NR<long>, FP_NR<dpe_t>>;
template class MatGSOInterface<Z_NR<double>, FP_NR<dpe_t>>;
template class MatGSOInterface<Z_NR<>, FP_NR<dpe_t>>;
template void adjust_radius_to_gh_bound<FP_NR<dpe_t>>(FP_NR<dpe_t> &max_dist, long max_dist_expo,
                                                      int block_size, const FP_NR<dpe_t> &root_det,
                                                      double gh_factor);
#endif

template class MatGSOInterface<Z_NR<long>, FP_NR<>>;
template class MatGSOInterface<Z_NR<double>, FP_NR<>>;
template class MatGSOInterface<Z_NR<>, FP_NR<>>;
template void adjust_radius_to_gh_bound<FP_NR<>>(FP_NR<> &max_dist, long max_dist_expo,
                                                       int block_size,
                                                       const FP_NR<> &root_det,
                                                       double gh_factor);

FPLLL_END_NAMESPACE
