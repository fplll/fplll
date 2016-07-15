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

#include "gso.h"

FPLLL_BEGIN_NAMESPACE

template <class ZT, class FT>
inline void MatGSO<ZT, FT>::invalidate_gso_row(int i, int new_valid_cols)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && new_valid_cols >= 0 && new_valid_cols <= i + 1);
  gso_valid_cols[i] = min(gso_valid_cols[i], new_valid_cols);
}

template <class ZT, class FT> void MatGSO<ZT, FT>::update_bf(int i)
{
  int n = max(n_known_cols, init_row_size[i]);
  if (enable_row_expo)
  {
    long max_expo = LONG_MIN;
    for (int j = 0; j < n; j++)
    {
      b(i, j).get_f_exp(bf(i, j), tmp_col_expo[j]);
      max_expo = max(max_expo, tmp_col_expo[j]);
    }
    for (int j = 0; j < n; j++)
    {
      bf(i, j).mul_2si(bf(i, j), tmp_col_expo[j] - max_expo);
    }
    row_expo[i] = max_expo;
  }
  else
  {
    for (int j = 0; j < n; j++)
    {
      bf(i, j).set_z(b(i, j));
    }
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::invalidate_gram_row(int i)
{
  for (int j = 0; j <= i; j++)
    gf(i, j).set_nan();
}

template <class ZT, class FT> void MatGSO<ZT, FT>::row_op_end(int first, int last)
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

template <class ZT, class FT> void MatGSO<ZT, FT>::discover_row()
{
  FPLLL_DEBUG_CHECK(n_known_rows < d);
  /* Early reduction (cols_locked=true) is not allowed when enable_int_gram=true,
     since n_known_cols might be too small to compute all the g(i,j). */
  FPLLL_DEBUG_CHECK(!(cols_locked && enable_int_gram));
  int i = n_known_rows;

  n_known_rows++;
  if (!cols_locked)
  {
    n_source_rows = n_known_rows;
    n_known_cols  = max(n_known_cols, init_row_size[i]);
  }
  if (enable_int_gram)
  {
    for (int j = 0; j <= i; j++)
      dot_product(g(i, j), b[i], b[j], n_known_cols);
  }
  else
  {
    invalidate_gram_row(i);
  }
  gso_valid_cols[i] = 0;
}

template <class ZT, class FT> inline ZT MatGSO<ZT, FT>::get_max_gram()
{
  ZT tmp;
  if (enable_int_gram)
  {
    tmp = g(0, 0);
    for (int i = 0; i < n_known_rows; i++)
      tmp = tmp.max_z(g(i, i));
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

template <class ZT, class FT> inline FT MatGSO<ZT, FT>::get_max_bstar()
{
  FT tmp;
  tmp = r(0, 0);
  for (int i = 0; i < n_known_rows; i++)
    tmp = tmp.max_f(r(i, i));
  return tmp;
}

template <class ZT, class FT> long MatGSO<ZT, FT>::get_max_mu_exp(int i, int n_columns)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && gso_valid_cols[i] >= n_columns);
  long max_expo = LONG_MIN, expo;
  for (int j = 0; j < n_columns; j++)
  {
    long expo2 = get_mu_exp(i, j, expo).exponent();
    max_expo    = max(max_expo, expo + expo2);
  }
  return max_expo;
}

template <class ZT, class FT> bool MatGSO<ZT, FT>::update_gso_row(int i, int last_j)
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

template <class ZT, class FT> void MatGSO<ZT, FT>::row_add(int i, int j)
{
  b[i].add(b[j], n_known_cols);
  if (enable_transform)
  {
    u[i].add(u[j]);
    if (enable_inverse_transform)
      u_inv_t[j].sub(u_inv_t[i]);
  }

  if (enable_int_gram)
  {
    // g(i, i) += 2 * g(i, j) + g(j, j)
    ztmp1.mul_2si(g(i, j), 1);
    ztmp1.add(ztmp1, g(j, j));
    g(i, i).add(g(i, i), ztmp1);

    for (int k = 0; k < n_known_rows; k++)
      if (k != i)
        sym_g(i, k).add(sym_g(i, k), sym_g(j, k));
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::row_sub(int i, int j)
{
  b[i].sub(b[j], n_known_cols);
  if (enable_transform)
  {
    u[i].sub(u[j]);
    if (enable_inverse_transform)
      u_inv_t[j].add(u_inv_t[i]);
  }

  if (enable_int_gram)
  {
    // g(i, i) += g(j, j) - 2 * g(i, j)
    ztmp1.mul_2si(g(i, j), 1);
    ztmp1.sub(g(j, j), ztmp1);
    g(i, i).add(g(i, i), ztmp1);

    for (int k = 0; k < n_known_rows; k++)
      if (k != i)
        sym_g(i, k).sub(sym_g(i, k), sym_g(j, k));
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::row_addmul_si(int i, int j, long x)
{
  b[i].addmul_si(b[j], x, n_known_cols);
  if (enable_transform)
  {
    u[i].addmul_si(u[j], x);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si(u_inv_t[i], -x);
  }

  if (enable_int_gram)
  {
    /* g(i, i) += 2 * (2^e * x) * g(i, j) + 2^(2*e) * x^2 * g(j, j)
      (must be done before updating g(i, j)) */
    ztmp1.mul_si(g(i, j), x);
    ztmp1.mul_2si(ztmp1, 1);
    g(i, i).add(g(i, i), ztmp1);
    ztmp1.mul_si(g(j, j), x);
    ztmp1.mul_si(ztmp1, x);
    g(i, i).add(g(i, i), ztmp1);

    // g(i, k) += g(j, k) * (2^e * x) for k != i
    for (int k = 0; k < n_known_rows; k++)
    {
      if (k == i)
        continue;
      ztmp1.mul_si(sym_g(j, k), x);
      sym_g(i, k).add(sym_g(i, k), ztmp1);
    }
  }
}

template <class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_si_2exp(int i, int j, long x, long expo)
{
  b[i].addmul_si_2exp(b[j], x, expo, n_known_cols, ztmp1);
  if (enable_transform)
  {
    u[i].addmul_si_2exp(u[j], x, expo, ztmp1);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si_2exp(u_inv_t[i], -x, expo, ztmp1);
  }

  if (enable_int_gram)
  {
    /* g(i, i) += 2 * (2^e * x) * g(i, j) + 2^(2*e) * x^2 * g(j, j)
      (must be done before updating g(i, j)) */
    ztmp1.mul_si(g(i, j), x);
    ztmp1.mul_2si(ztmp1, expo + 1);
    g(i, i).add(g(i, i), ztmp1);
    ztmp1.mul_si(g(j, j), x);
    ztmp1.mul_si(ztmp1, x);
    ztmp1.mul_2si(ztmp1, 2 * expo);
    g(i, i).add(g(i, i), ztmp1);

    // g(i, k) += g(j, k) * (2^e * x) for k != i
    for (int k = 0; k < n_known_rows; k++)
    {
      if (k == i)
        continue;
      ztmp1.mul_si(sym_g(j, k), x);
      ztmp1.mul_2si(ztmp1, expo);
      sym_g(i, k).add(sym_g(i, k), ztmp1);
    }
  }
}

template <class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_2exp(int i, int j, const ZT &x, long expo)
{
  b[i].addmul_2exp(b[j], x, expo, n_known_cols, ztmp1);
  if (enable_transform)
  {
    u[i].addmul_2exp(u[j], x, expo, ztmp1);
    if (enable_inverse_transform)
    {
      ZT minusX;
      minusX.neg(x);
      u_inv_t[j].addmul_2exp(u_inv_t[i], minusX, expo, ztmp1);
    }
  }

  if (enable_int_gram)
  {
    /* g(i, i) += 2 * (2^e * x) * g(i, j) + 2^(2*e) * x^2 * g(j, j)
      (must be done before updating g(i, j)) */
    ztmp1.mul(g(i, j), x);
    ztmp1.mul_2si(ztmp1, expo + 1);
    g(i, i).add(g(i, i), ztmp1);
    ztmp1.mul(g(j, j), x);
    ztmp1.mul(ztmp1, x);
    ztmp1.mul_2si(ztmp1, 2 * expo);
    g(i, i).add(g(i, i), ztmp1);

    // g(i, k) += g(j, k) * (2^e * x) for k != i
    for (int k = 0; k < n_known_rows; k++)
    {
      if (k == i)
        continue;
      ztmp1.mul(sym_g(j, k), x);
      ztmp1.mul_2si(ztmp1, expo);
      sym_g(i, k).add(sym_g(i, k), ztmp1);
    }
  }
}

template <class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_we(int i, int j, const FT &x, long expo_add)
{
  FPLLL_DEBUG_CHECK(j >= 0 && /*i > j &&*/ i < n_known_rows && j < n_source_rows);
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

template <class ZT, class FT> void MatGSO<ZT, FT>::row_swap(int i, int j)
{
  FPLLL_DEBUG_CHECK(!enable_inverse_transform);
  b.swapRows(i, j);
  if (enable_transform)
  {
    u.swapRows(i, j);
  }

  if (enable_int_gram)
  {
    for (int k = 0; k < i; k++)
      g(i, k).swap(g(j, k));
    for (int k = i + 1; k < j; k++)
      g(k, i).swap(g(j, k));
    for (int k = j + 1; k < n_known_rows; k++)
      g(k, i).swap(g(k, j));
    g(i, i).swap(g(j, j));
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::move_row(int old_r, int new_r)
{
  FPLLL_DEBUG_CHECK(!cols_locked);
  if (new_r < old_r)
  {
    FPLLL_DEBUG_CHECK(old_r < n_known_rows && !cols_locked);
    for (int i = new_r; i < n_known_rows; i++)
    {
      invalidate_gso_row(i, new_r);
    }
    rotate(gso_valid_cols.begin() + new_r, gso_valid_cols.begin() + old_r,
           gso_valid_cols.begin() + old_r + 1);
    mu.rotate_right(new_r, old_r);
    r.rotate_right(new_r, old_r);
    b.rotate_right(new_r, old_r);
    if (enable_transform)
    {
      u.rotate_right(new_r, old_r);
      if (enable_inverse_transform)
        u_inv_t.rotate_right(new_r, old_r);
    }
    if (enable_int_gram)
      g.rotate_gram_right(new_r, old_r, n_known_rows);
    else
    {
      gf.rotate_gram_right(new_r, old_r, n_known_rows);
      bf.rotate_right(new_r, old_r);
    }
    if (enable_row_expo)
      rotate(row_expo.begin() + new_r, row_expo.begin() + old_r, row_expo.begin() + old_r + 1);
  }
  else if (new_r > old_r)
  {
    for (int i = old_r; i < n_known_rows; i++)
    {
      invalidate_gso_row(i, old_r);
    }
    rotate(gso_valid_cols.begin() + old_r, gso_valid_cols.begin() + old_r + 1,
           gso_valid_cols.begin() + new_r + 1);
    mu.rotate_left(old_r, new_r);
    r.rotate_left(old_r, new_r);
    b.rotate_left(old_r, new_r);
    if (enable_transform)
    {
      u.rotate_left(old_r, new_r);
      if (enable_inverse_transform)
        u_inv_t.rotate_left(old_r, new_r);
    }
    if (enable_int_gram)
    {
      if (old_r < n_known_rows - 1)
        g.rotate_gram_left(old_r, min(new_r, n_known_rows - 1), n_known_rows);
    }
    else
    {
      if (old_r < n_known_rows - 1)
        gf.rotate_gram_left(old_r, min(new_r, n_known_rows - 1), n_known_rows);
      bf.rotate_left(old_r, new_r);
    }
    if (enable_row_expo)
      rotate(row_expo.begin() + old_r, row_expo.begin() + old_r + 1, row_expo.begin() + new_r + 1);
    if (new_r >= n_known_rows)
    {
      rotate(init_row_size.begin() + old_r, init_row_size.begin() + old_r + 1,
             init_row_size.begin() + new_r + 1);
      if (old_r < n_known_rows)
      {
        n_known_rows--;
        n_source_rows        = n_known_rows;
        init_row_size[new_r] = max(b[new_r].size_nz(), 1);
      }
    }
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::lock_cols() { cols_locked = true; }

template <class ZT, class FT> void MatGSO<ZT, FT>::unlock_cols()
{
  n_known_rows = n_source_rows;
  cols_locked  = false;
}

template <class ZT, class FT>
void MatGSO<ZT, FT>::apply_transform(const Matrix<FT> &transform, int src_base, int target_base)
{
  int targetSize = transform.get_rows(), srcSize = transform.get_cols();
  int oldD = d;
  create_rows(targetSize);
  for (int i = 0; i < targetSize; i++)
  {
    for (int j = 0; j < srcSize; j++)
    {
      row_addmul(oldD + i, src_base + j, transform(i, j));
    }
  }
  row_op_begin(target_base, target_base + targetSize);
  for (int i = 0; i < targetSize; i++)
  {
    row_swap(target_base + i, oldD + i);
  }
  row_op_end(target_base, target_base + targetSize);
  remove_last_rows(targetSize);
}

template <class ZT, class FT> void MatGSO<ZT, FT>::size_increased()
{
  int oldD = mu.get_rows();

  if (d > alloc_dim)
  {
    if (enable_int_gram)
      g.resize(d, d);
    else
    {
      bf.resize(d, b.get_cols());
      gf.resize(d, d);
    }
    mu.resize(d, d);
    r.resize(d, d);
    gso_valid_cols.resize(d);
    init_row_size.resize(d);
    if (enable_row_expo)
    {
      row_expo.resize(d);
    }
    alloc_dim = d;
  }

  for (int i = oldD; i < d; i++)
  {
    init_row_size[i] = max(b[i].size_nz(), 1);
    if (!enable_int_gram)
    {
      bf[i].fill(0);  // update_bf might not copy all the zeros of b[i]
      update_bf(i);
    }
  }
}

template <class ZT, class FT> double MatGSO<ZT, FT>::get_current_slope(int start_row, int stop_row)
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

template <class ZT, class FT> FT MatGSO<ZT, FT>::get_root_det(int start_row, int end_row)
{
  start_row   = max(0, start_row);
  end_row     = min(d, end_row);
  FT h        = (double)(end_row - start_row);
  FT root_det = get_log_det(start_row, end_row) / h;
  root_det.exponential(root_det);
  return root_det;
}

template <class ZT, class FT> FT MatGSO<ZT, FT>::get_log_det(int start_row, int end_row)
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
FT MatGSO<ZT, FT>::get_slide_potential(int start_row, int end_row, int block_size)
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
void gaussian_heuristic(FT &max_dist, long max_dist_expo, int block_size, const FT &root_det,
                        double gh_factor)
{
  double t = (double)block_size / 2.0 + 1;
  t        = tgamma(t);
  t        = pow(t, 2.0 / (double)block_size);
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

template class MatGSO<Z_NR<long>, FP_NR<double>>;
template class MatGSO<Z_NR<double>, FP_NR<double>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<double>>;
template void gaussian_heuristic<FP_NR<double>>(FP_NR<double> &max_dist, long max_dist_expo,
                                                int block_size, const FP_NR<double> &root_det,
                                                double gh_factor);

#ifdef FPLLL_WITH_LONG_DOUBLE
template class MatGSO<Z_NR<long>, FP_NR<long double>>;
template class MatGSO<Z_NR<double>, FP_NR<long double>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<long double>>;
template void gaussian_heuristic<FP_NR<long double>>(FP_NR<long double> &max_dist,
                                                     long max_dist_expo, int block_size,
                                                     const FP_NR<long double> &root_det,
                                                     double gh_factor);

#endif

#ifdef FPLLL_WITH_QD
template class MatGSO<Z_NR<long>, FP_NR<dd_real>>;
template class MatGSO<Z_NR<double>, FP_NR<dd_real>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<dd_real>>;
template void gaussian_heuristic<FP_NR<dd_real>>(FP_NR<dd_real> &max_dist, long max_dist_expo,
                                                 int block_size, const FP_NR<dd_real> &root_det,
                                                 double gh_factor);
template class MatGSO<Z_NR<long>, FP_NR<qd_real>>;
template class MatGSO<Z_NR<double>, FP_NR<qd_real>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<qd_real>>;
template void gaussian_heuristic<FP_NR<qd_real>>(FP_NR<qd_real> &max_dist, long max_dist_expo,
                                                 int block_size, const FP_NR<qd_real> &root_det,
                                                 double gh_factor);
#endif

#ifdef FPLLL_WITH_DPE
template class MatGSO<Z_NR<long>, FP_NR<dpe_t>>;
template class MatGSO<Z_NR<double>, FP_NR<dpe_t>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<dpe_t>>;
template void gaussian_heuristic<FP_NR<dpe_t>>(FP_NR<dpe_t> &max_dist, long max_dist_expo,
                                               int block_size, const FP_NR<dpe_t> &root_det,
                                               double gh_factor);
#endif

template class MatGSO<Z_NR<long>, FP_NR<mpfr_t>>;
template class MatGSO<Z_NR<double>, FP_NR<mpfr_t>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>>;
template void gaussian_heuristic<FP_NR<mpfr_t>>(FP_NR<mpfr_t> &max_dist, long max_dist_expo,
                                                int block_size, const FP_NR<mpfr_t> &root_det,
                                                double gh_factor);

FPLLL_END_NAMESPACE
