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

#ifndef FPLLL_GSO_HOUSEHOLDER_H
#define FPLLL_GSO_HOUSEHOLDER_H

#include "gso_interface.h"
#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

/**
 * MatGSOHouseholder provides an interface for performing elementary operations on a basis
 * and computing its Gram matrix and its Gram-Schmidt orthogonalization.
 * The Gram-Schmidt coefficients are computed on demand. The object keeps track
 * of which coefficients are valid after each row operation.
 */
template <class ZT, class FT> class MatGSOHouseholder : public MatGSOInterface<ZT, FT>
{
public:
  using MatGSOInterface<ZT, FT>::d;
  using MatGSOInterface<ZT, FT>::n_known_rows;
  using MatGSOInterface<ZT, FT>::n_source_rows;
  using MatGSOInterface<ZT, FT>::cols_locked;  // maybe scratch.
  using MatGSOInterface<ZT, FT>::enable_int_gram;
  using MatGSOInterface<ZT, FT>::gso_valid_cols;
  using MatGSOInterface<ZT, FT>::alloc_dim;
  using MatGSOInterface<ZT, FT>::bf;
  using MatGSOInterface<ZT, FT>::ftmp1;
  using MatGSOInterface<ZT, FT>::ftmp2;
  using MatGSOInterface<ZT, FT>::init_row_size;
  using MatGSOInterface<ZT, FT>::enable_row_expo;
  using MatGSOInterface<ZT, FT>::row_expo;
  using MatGSOInterface<ZT, FT>::n_known_cols;
  using MatGSOInterface<ZT, FT>::tmp_col_expo;
  using MatGSOInterface<ZT, FT>::update_gso;
  using MatGSOInterface<ZT, FT>::update_gso_row;

#ifdef DEBUG
  /* Used only in debug mode. */
  using MatGSOInterface<ZT, FT>::row_op_first;
  using MatGSOInterface<ZT, FT>::row_op_last;
  using MatGSOInterface<ZT, FT>::in_row_op_range;
#endif

  /**
   * Constructor.
   * The precision of FT must be defined before creating an instance of the
   * class and must remain the same until the object is destroyed (or no longer
   * needed).
   * @param b
   *   The matrix on which row operations are performed. It must not be empty.
   * @param u
   *   If u is not empty, operations on b are also done on u
   *   (in this case both must have the same number of rows).
   *   If u is initially the identity matrix, multiplying transform by the
   *   initial basis gives the current basis.
   * @param u_inv_t
   *   Inverse transform (should be empty, which disables the computation, or
   *   initialized with identity matrix). It works only if u is not empty.
   * @param enable_int_gram
   *   If true, coefficients of the Gram matrix are computed with exact integer
   *   arithmetic (type ZT). Otherwise, they are computed in floating-point
   *   (type FT). Note that when exact arithmetic is used, all coefficients of
   *   the first n_known_rows are continuously updated, whereas in floating-point,
   *   they are computed only on-demand. This option cannot be enabled if
   *   enable_row_expo=true.
   * @param enable_row_expo
   *   If true, each row of b is normalized by a power of 2 before doing
   *   conversion to floating-point, which hopefully avoids some overflows.
   *   This option cannot be enabled if enable_int_gram=true and works only
   *   with FT=double and FT=long double. It is useless and MUST NOT be used
   *   for FT=dpe or FT=mpfr_t.
   * @param row_op_force_long
   *   Affects the behaviour of row_addmul(_we).
   *   See the documentation of row_addmul.
   */
  MatGSOHouseholder(Matrix<ZT> &arg_b, Matrix<ZT> &arg_u, Matrix<ZT> &arg_uinv_t, int flags)
      : MatGSOInterface<ZT, FT>(arg_u, arg_uinv_t, flags), b(arg_b)
  {
    FPLLL_DEBUG_CHECK(!(enable_int_gram && enable_row_expo));
    d = b.get_rows();
    n = b.get_cols();
    if (enable_row_expo)
    {
      tmp_col_expo.resize(b.get_cols());
    }
    size_increased();
#ifdef DEBUG
    row_op_first = row_op_last = -1;
#endif
  }

  ~MatGSOHouseholder() { delete[] sigma_householder; }

public:
  /**
   * Basis of the lattice
   */
  Matrix<ZT> &b;

  virtual inline long get_max_exp_of_b();
  virtual inline bool b_row_is_zero(int i);
  virtual inline int get_cols_of_b();
  virtual inline int get_rows_of_b();
  virtual inline void negate_row_of_b(int i);

  //  virtual inline void printparam(ostream &os);
  virtual inline FT &get_gram(FT &f, int i, int j);

  virtual void move_row(int old_r, int new_r) {}
  virtual void row_addmul_we(int i, int j, const FT &x, long expo_add) {}
  virtual void row_add(int i, int j) {}
  virtual void row_sub(int i, int j) {}
  virtual void create_rows(int n_new_rows) {}
  virtual void remove_last_rows(int n_removed_rows) {}
  virtual void row_swap(int i, int j) {}
  virtual bool update_gso_row(int i, int last_j);

  /**
   * Returns f = householder_r(i, j).
   *
   * Returns reference to `f`.
   */
  inline FT &get_r_householder(FT &f, int i, int j);

  /**
   * Returns the r_householder matrix
   */
  const Matrix<FT> &get_r_householder_matrix() { return r_householder; }

  /**
   * Apply Householder transformation on row i.
   */
  void update_r_householder_row(int i);

private:
  /* Allocates matrices and arrays whose size depends on d (all but tmp_col_expo).
   When enable_int_gram=false, initializes bf. */
  virtual void size_increased();

  virtual void discover_row();

  /* Upates the i-th row of bf. It does not invalidate anything, so the caller
     must take into account that it might change row_expo. */
  virtual void update_bf(int i);
  virtual void invalidate_gram_row(int i) {}
  virtual void row_addmul_si(int i, int j, long x) {}

  // Number of columns of b.
  int n;

  /**
   * bf = r_householder * q_householder.
   * r_householder is lower triangular and the diagonal coefficient are >= 0.
   */
  Matrix<FT> r_householder;

  /**
   * Vector v following [MSV, ISSAC'09].
   */
  Matrix<FT> v_householder;

  /**
   * Sigma values following [MSV, ISSAC'09].
   */
  FT *sigma_householder;
};

template <class ZT, class FT> inline long MatGSOHouseholder<ZT, FT>::get_max_exp_of_b()
{
  return b.get_max_exp();
}

template <class ZT, class FT> inline bool MatGSOHouseholder<ZT, FT>::b_row_is_zero(int i)
{
  return b[i].is_zero();
}
template <class ZT, class FT> inline int MatGSOHouseholder<ZT, FT>::get_cols_of_b()
{
  return b.get_cols();
}

template <class ZT, class FT> inline int MatGSOHouseholder<ZT, FT>::get_rows_of_b()
{
  return b.get_rows();
}

template <class ZT, class FT> inline void MatGSOHouseholder<ZT, FT>::negate_row_of_b(int i)
{

  for (int j = 0; j < get_cols_of_b(); j++)
  {
    b[i][j].neg(b[i][j]);
  }
}

template <class ZT, class FT> inline FT &MatGSOHouseholder<ZT, FT>::get_gram(FT &f, int i, int j)
{
  f.set_nan();
  return f;
}

template <class ZT, class FT>
inline FT &MatGSOHouseholder<ZT, FT>::get_r_householder(FT &f, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && j >= 0 && j <= i);
  f = r_householder(i, j);
  return f;
}

FPLLL_END_NAMESPACE

#endif  // FPLLL_GSO_HOUSEHOLDER_H
