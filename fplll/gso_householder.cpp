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

#include "gso_householder.h"

FPLLL_BEGIN_NAMESPACE

template <class ZT, class FT> void MatGSOHouseholder<ZT, FT>::update_bf(int i)
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

template <class ZT, class FT> void MatGSOHouseholder<ZT, FT>::discover_row()
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

  update_r_householder_row(i);

  gso_valid_cols[i] = 0;
}

template <class ZT, class FT> bool MatGSOHouseholder<ZT, FT>::update_gso_row(int i, int last_j)
{
  // FPLLL_TRACE_IN("Updating GSO up to (" << i << ", " << last_j << ")");
  // FPLLL_TRACE("n_known_rows=" << n_known_rows << " n_source_rows=" << n_source_rows);
  if (i >= n_known_rows)
  {
    discover_row();
  }
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && last_j >= 0 && last_j < n_source_rows);

  int j = max(0, gso_valid_cols[i]);

  gso_valid_cols[i] = j;  // = max(0, gso_valid_cols[i], last_j + 1)
  // FPLLL_TRACE_OUT("End of GSO update");
  return true;
}

#if 0
template <class ZT, class FT> void MatGSOHouseholder<ZT, FT>::update_r_householder_row(int i)
{
  FT norm_r_square;
  dot_product(norm_r_square, r_householder[i], r_householder[i], i, n);

  if (r_householder(i, i) * r_householder(i, i) != norm_r_square)
  {
    FT snorm_r = sqrt(norm_r_square);
    if (r_householder(i, i).cmp(0) > 0)
    {
      snorm_r = -snorm_r;
    }
    // u = r + snorm_r * e_i
    FT norm_u_square = norm_r_square + norm_r_square + 2.0 * r_householder(i, i) * snorm_r;
    ftmp1            = 2.0 / norm_u_square;

    for (int j = i + 1; j < d; j++)
    {
      dot_product(ftmp2, r_householder[i], r_householder[j], i, n);
      ftmp2.addmul(r_householder(j, i), snorm_r);
      // At this point, dot_product(ftmp2, u, r_householder[j], i, n)
      ftmp2 *= ftmp1;

      // r_householder(j, i) = r_householder[j][i] - dot_product(u, R[j], i, c) * ftmp1 * u[k]
      r_householder(j, i).submul(ftmp2, r_householder(i, i) + snorm_r);

      for (int k = i + 1; k < n; k++)
      {
        // r_householder(j, k) = r_householder(j, k) - dot_product(u, R[j], i, c) * ftmp1 * u[k]
        r_householder(j, k).submul(ftmp2, r_householder(i, k));
      }
    }

    ftmp2 = norm_r_square;
    ftmp2.addmul(r_householder(i, i), snorm_r);
    // At this point, dot_product(ftmp2, u, r_householder[i], i, n
    ftmp2 *= ftmp1;
    // r_householder(i, i) = r_householder(i, i) - dot_product(u, r_householder[i], i, c) * ftmp1
    // * u[i]
    r_householder(i, i).submul(ftmp2, r_householder(i, i) + snorm_r);

    for (int k = i + 1; k < n; k++)
    {
      r_householder(i, k) = 0.0;
      FPLLL_DEBUG_CHECK(r_householder(i, k).is_zero());
    }
  }

  // r_householder(i, i) must be non-negative
  if (r_householder(i, i).cmp(0) < 0)
  {
    for (int k = i; k < d; k++)
    {
      r_householder(k, i) = -r_householder(k, i);
    }
  }
}
#else   // 0
template <class ZT, class FT> void MatGSOHouseholder<ZT, FT>::update_r_householder_row(int i)
{
  int j, k;
  if (!enable_int_gram)
  {
    for (j = 0; j < n; j++)
    {
      r_householder(i, j) = bf(i, j);
    }
  }
  else
  {
    for (j = 0; j < n; j++)
    {
      r_householder(i, j).set_z(b(i, j));
    }
  }
  FT tmp;
  for (j = 0; j < i; j++)
  {
    // vj * ri[j..n]^T
    dot_product(ftmp1, v_householder[j], r_householder[i], j, n);
    //-vj * ri[j..n]^T
    ftmp1.neg(ftmp1);
    for (k = j; k < n; k++)
    {
      // ri[j..n] = ri[j..n] - (vj * ri[j..n]^T) * vj
      r_householder(i, k).addmul(v_householder(j, k), ftmp1);
    }
    // ri[j] = sigma[j] * ri[j]
    r_householder(i, j).mul(sigma_householder[j], r_householder(i, j));
  }
  // sigma[i] = sign(r[1])
  sigma_householder[i] = (r_householder(i, i).cmp(0) < 0) ? -1.0 : 1.0;
  // r^T * r
  dot_product(ftmp1, r_householder[i], r_householder[i], i, n);
  if (ftmp1.cmp(0) != 0)
  {
    ftmp2.sqrt(ftmp1);
    // s = sigma[i] * ||r|| = sigma[i] * sqrt(r * r^T)
    tmp.mul(sigma_householder[i], ftmp2);
    v_householder(i, i).mul(r_householder(i, i), r_householder(i, i));
    v_householder(i, i).sub(v_householder(i, i), ftmp1);
    ftmp1.add(r_householder(i, i), tmp);
    v_householder(i, i).div(v_householder(i, i), ftmp1);
    // Here, vi[1] = (-sum(r[j]^2, j, 2, n-i+1) / (r[1] + s)
    if (v_householder(i, i).cmp(0) != 0)
    {
      tmp.neg(tmp);
      tmp.mul(tmp, v_householder(i, i));
      tmp.sqrt(tmp);
      tmp.div(1.0, tmp);
      // Here, tmp = 1 / sqrt(-s * vi[1])
      v_householder(i, i).mul(v_householder(i, i), tmp);
      r_householder(i, i) = ftmp2;
      for (k = i + 1; k < n; k++)
      {
        v_householder(i, k).mul(r_householder(i, k), tmp);
        r_householder(i, k) = 0.0;
        FPLLL_DEBUG_CHECK(r_householder(i, k).is_zero());
      }
      // Here, vi = vi / tmp and ri[i..n] = (||r||, 0, 0, ..., 0)
    }
  }
}
#endif  // 0

template <class ZT, class FT> void MatGSOHouseholder<ZT, FT>::size_increased()
{
  int old_d = bf.get_rows();

  if (d > alloc_dim)
  {
    if (!enable_int_gram)
    {
      bf.resize(d, b.get_cols());
    }
    gso_valid_cols.resize(d);
    init_row_size.resize(d);
    if (enable_row_expo)
    {
      row_expo.resize(d);
    }
    sigma_householder = new FT[b.get_rows()];
    r_householder.resize(d, n);
    v_householder.resize(d, n);
#ifdef DEBUG
    for (int i = 0; i < d; i++)
    {
      for (int j = 0; j < n; j++)
      {
        v_householder(i, j).set_nan();
      }
    }
#endif  // DEBUG
    alloc_dim = d;
  }

  for (int i = old_d; i < d; i++)
  {
    init_row_size[i] = max(b[i].size_nz(), 1);
    if (!enable_int_gram)
    {
      bf[i].fill(0);  // update_bf might not copy all the zeros of b[i]
      update_bf(i);
    }
  }
}

template class MatGSOHouseholder<Z_NR<long>, FP_NR<double>>;
template class MatGSOHouseholder<Z_NR<double>, FP_NR<double>>;
template class MatGSOHouseholder<Z_NR<mpz_t>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class MatGSOHouseholder<Z_NR<long>, FP_NR<long double>>;
template class MatGSOHouseholder<Z_NR<double>, FP_NR<long double>>;
template class MatGSOHouseholder<Z_NR<mpz_t>, FP_NR<long double>>;

#endif

#ifdef FPLLL_WITH_QD
template class MatGSOHouseholder<Z_NR<long>, FP_NR<dd_real>>;
template class MatGSOHouseholder<Z_NR<double>, FP_NR<dd_real>>;
template class MatGSOHouseholder<Z_NR<mpz_t>, FP_NR<dd_real>>;

template class MatGSOHouseholder<Z_NR<long>, FP_NR<qd_real>>;
template class MatGSOHouseholder<Z_NR<double>, FP_NR<qd_real>>;
template class MatGSOHouseholder<Z_NR<mpz_t>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class MatGSOHouseholder<Z_NR<long>, FP_NR<dpe_t>>;
template class MatGSOHouseholder<Z_NR<double>, FP_NR<dpe_t>>;
template class MatGSOHouseholder<Z_NR<mpz_t>, FP_NR<dpe_t>>;
#endif

template class MatGSOHouseholder<Z_NR<long>, FP_NR<mpfr_t>>;
template class MatGSOHouseholder<Z_NR<double>, FP_NR<mpfr_t>>;
template class MatGSOHouseholder<Z_NR<mpz_t>, FP_NR<mpfr_t>>;

FPLLL_END_NAMESPACE
