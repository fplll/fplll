/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2011 Xavier Pujol.
   Copyright (C) 2019 Koen de Boer

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

#ifndef FPLLL_GSOGRAM_H
#define FPLLL_GSOGRAM_H

#include "gso_interface.h"
#include "nr/matrix.h"
#include "util.h"

FPLLL_BEGIN_NAMESPACE

/**
 * MatGSOGram provides an interface for performing elementary operations on a basis
 * and computing its Gram matrix and its Gram-Schmidt orthogonalization.
 * The Gram-Schmidt coefficients are computed on demand. The object keeps track
 * of which coefficients are valid after each row operation.
 */
template <class ZT, class FT> class MatGSOGram : public MatGSOInterface<ZT, FT>
{
public:
  using MatGSOInterface<ZT, FT>::d;
  using MatGSOInterface<ZT, FT>::n_known_rows;
  using MatGSOInterface<ZT, FT>::n_source_rows;
  using MatGSOInterface<ZT, FT>::u;
  using MatGSOInterface<ZT, FT>::enable_transform;
  using MatGSOInterface<ZT, FT>::cols_locked;  // maybe scratch.
  using MatGSOInterface<ZT, FT>::enable_int_gram;
  using MatGSOInterface<ZT, FT>::gso_valid_cols;
  using MatGSOInterface<ZT, FT>::enable_inverse_transform;
  using MatGSOInterface<ZT, FT>::u_inv_t;
  using MatGSOInterface<ZT, FT>::sym_g;
  using MatGSOInterface<ZT, FT>::mu;
  using MatGSOInterface<ZT, FT>::r;
  using MatGSOInterface<ZT, FT>::ztmp1;
  using MatGSOInterface<ZT, FT>::ztmp2;
  using MatGSOInterface<ZT, FT>::row_op_force_long;
  using MatGSOInterface<ZT, FT>::alloc_dim;
  using MatGSOInterface<ZT, FT>::get_mu;
  using MatGSOInterface<ZT, FT>::get_r;
  using MatGSOInterface<ZT, FT>::gptr;
  using MatGSOInterface<ZT, FT>::invalidate_gso_row;

  using MatGSOInterface<ZT, FT>::create_row;
  using MatGSOInterface<ZT, FT>::remove_last_row;
  using MatGSOInterface<ZT, FT>::print_mu_r_g;
  using MatGSOInterface<ZT, FT>::discover_all_rows;
  using MatGSOInterface<ZT, FT>::update_gso;
  using MatGSOInterface<ZT, FT>::update_gso_row;
  using MatGSOInterface<ZT, FT>::symmetrize_g;

#ifdef DEBUG
  /* Used only in debug mode. */
  using MatGSOInterface<ZT, FT>::row_op_first;
  using MatGSOInterface<ZT, FT>::row_op_last;
  using MatGSOInterface<ZT, FT>::in_row_op_range;
#endif

  MatGSOGram(Matrix<ZT> &arg_g, Matrix<ZT> &arg_u, Matrix<ZT> &arg_uinv_t, int flags = GSO_INT_GRAM)
      : MatGSOInterface<ZT, FT>(arg_u, arg_uinv_t, flags)
  {
    if (flags != GSO_INT_GRAM)
    {
      throw std::invalid_argument("flags must be equal to GSO_INT_GRAM");
    }
    gptr = &arg_g;
    if (gptr == nullptr)
    {
      throw std::runtime_error("Error: gptr is equal to the nullpointer.");
    }
    d = gptr->get_rows();
    size_increased();

#ifdef DEBUG
    row_op_first = row_op_last = -1;
#endif
  }

public:
  virtual inline long get_max_exp_of_b();
  virtual inline bool b_row_is_zero(int i);
  virtual inline int get_cols_of_b() const;
  virtual inline int get_rows_of_b() const;
  virtual inline void negate_row_of_b(int i);

  //  inline void set_g(Matrix<ZT> arg_g)
  virtual inline void create_rows(int n_new_rows);
  virtual inline void remove_last_rows(int n_removed_rows);

  /*
    virtual inline void dump_mu_d(double *mu, int offset = 0, int block_size = -1);
    virtual inline void dump_mu_d(vector<double> mu, int offset = 0, int block_size = -1);

    virtual inline void dump_r_d(double *r, int offset = 0, int block_size = -1);
    virtual inline void dump_r_d(vector<double> r, int offset = 0, int block_size = -1);
  */
  virtual void move_row(int old_r, int new_r);

  virtual void row_addmul_we(int i, int j, const FT &x, long expo_add);

  // b[i] += b[j] / b[i] -= b[j] (i > j)
  virtual void row_add(int i, int j);
  virtual void row_sub(int i, int j);

  //  virtual inline void printparam(ostream &os);

  virtual inline ZT &sqnorm_coordinates(ZT &sqnorm, vector<ZT> coordinates);

  virtual inline FT &get_gram(FT &f, int i, int j);

  virtual inline ZT &get_int_gram(ZT &z, int i, int j);

  // b[i] <-> b[j] (i < j)
  virtual void row_swap(int i, int j);

private:
  /* Allocates matrices and arrays whose size depends on d (all but tmp_col_expo).
   When enable_int_gram=false, initializes bf. */
  virtual void size_increased();

  virtual void discover_row();

  /* Upates the i-th row of bf. It does not invalidate anything, so the caller
     must take into account that it might change row_expo. */
  virtual void update_bf(int i);
  /* Marks g(i, j) for all j <= i (but NOT for j > i) */
  virtual void invalidate_gram_row(int i);

  // b[i] <- b[i] + x * b[j] (i > j)
  virtual void row_addmul_si(int i, int j, long x);
  // b[i] <- b[i] + (2^expo * x) * b[j] (i > j)

  virtual void row_addmul_si_2exp(int i, int j, long x, long expo);
  virtual void row_addmul_2exp(int i, int j, const ZT &x, long expo);
  // virtual void apply_transform(const Matrix<FT> &transform, int src_base, int target_base);
};

template <class ZT, class FT>
inline ZT &MatGSOGram<ZT, FT>::sqnorm_coordinates(ZT &sqnorm, vector<ZT> coordinates)
{
  vector<ZT> tmpvec;
  Matrix<ZT> &g = *gptr;
  vector_matrix_product(tmpvec, coordinates, g);

  sqnorm = 0;
  for (int i = 0; i < g.get_cols(); i++)
  {
    ztmp1.mul(tmpvec[i], coordinates[i]);
    sqnorm.add(sqnorm, ztmp1);
  }
  return sqnorm;
}

template <class ZT, class FT> inline long MatGSOGram<ZT, FT>::get_max_exp_of_b()
{
  if (gptr == nullptr)
  {
    throw std::runtime_error("Error: gptr is equal to the nullpointer.");
  }
  Matrix<ZT> &g = *gptr;
  // normally, this returns the maximum
  // exponent of b. We approximate this temporarily
  // by the half of the maximum exponent of g.
  // beware of errors, because maybe g is
  // zero in the upper half of the matrix.
  return g.get_max_exp() / 2;
}
template <class ZT, class FT> inline bool MatGSOGram<ZT, FT>::b_row_is_zero(int i)
{
  if (gptr == nullptr)
  {
    throw std::runtime_error("Error: gptr is equal to the nullpointer.");
  }
  Matrix<ZT> &g = *gptr;
  // normally this returns whether the
  // ith row of the basis is the zero row.
  // now, we just check whether g[i][i] is zero.
  return g[i][i].is_zero();
}
template <class ZT, class FT> inline int MatGSOGram<ZT, FT>::get_cols_of_b() const
{
  if (gptr == nullptr)
  {
    throw std::runtime_error("Error: gptr is equal to the nullpointer.");
  }
  // in this gram-matrix version, the number
  // of columns of b is the same as the
  // number of colums of g.
  return gptr->get_cols();
}

template <class ZT, class FT> inline int MatGSOGram<ZT, FT>::get_rows_of_b() const
{
  if (gptr == nullptr)
  {
    throw std::runtime_error("Error: gptr is equal to the nullpointer.");
  }
  // in this gram-matrix version, the number
  // of columns of b is the same as the
  // number of colums of g.
  return gptr->get_rows();
}

template <class ZT, class FT> inline void MatGSOGram<ZT, FT>::negate_row_of_b(int i)
{
  if (enable_int_gram)
  {

    for (int j = 0; j < get_rows_of_b(); j++)
    {
      if (j != i)
      {
        sym_g(i, j).neg(sym_g(i, j));
      }
    }
  }
}

template <class ZT, class FT> inline FT &MatGSOGram<ZT, FT>::get_gram(FT &f, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && j >= 0 && j <= i && j < n_source_rows &&
                    !in_row_op_range(i));
  if (enable_int_gram)
  {
    if (gptr == nullptr)
    {
      throw std::runtime_error("Error: gptr is equal to the nullpointer.");
    }
    f.set_z((*gptr)(i, j));
  }
  return f;
}

template <class ZT, class FT> inline ZT &MatGSOGram<ZT, FT>::get_int_gram(ZT &z, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && j >= 0 && j <= i && j < n_source_rows &&
                    !in_row_op_range(i));
  if (enable_int_gram)
  {
    if (gptr == nullptr)
    {
      throw std::runtime_error("Error: gptr is equal to the nullpointer.");
    }
    z = (*gptr)[i][j];
  }
  return z;
}

template <class ZT, class FT> inline void MatGSOGram<ZT, FT>::remove_last_rows(int n_removed_rows)
{
  FPLLL_DEBUG_CHECK(!cols_locked && d >= n_removed_rows);
  d -= n_removed_rows;
  n_known_rows  = min(n_known_rows, d);
  n_source_rows = n_known_rows;
  if (enable_transform)
    u.set_rows(d);
}

template <class ZT, class FT> inline void MatGSOGram<ZT, FT>::create_rows(int n_new_rows)
{
  FPLLL_DEBUG_CHECK(!cols_locked);
  int old_d = d;
  d += n_new_rows;

  if (enable_transform)
  {
    u.set_rows(d);
    for (int i = old_d; i < d; i++)
      for (int j = 0; j < u.get_cols(); j++)
        u[i][j] = 0;
  }
  size_increased();
  if (n_known_rows == old_d)
    MatGSOInterface<ZT, FT>::discover_all_rows();
}

FPLLL_END_NAMESPACE

#endif
