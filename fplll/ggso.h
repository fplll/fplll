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

#ifndef FPLLL_GGSO_H
#define FPLLL_GGSO_H


#include "nr/matrix.h"
#include "gsob.h"

FPLLL_BEGIN_NAMESPACE

/**
 * MatGGSO provides an interface for performing elementary operations on a basis
 * and computing its Gram matrix and its Gram-Schmidt orthogonalization.
 * The Gram-Schmidt coefficients are computed on demand. The object keeps track
 * of which coefficients are valid after each row operation.
 */
template <class ZT, class FT> class MatGGSO : public MatGSOb<ZT,FT>
{
public:
  
  using MatGSOb<ZT,FT>::d;
  using MatGSOb<ZT,FT>::n_known_rows;
  using MatGSOb<ZT,FT>::n_source_rows;
  using MatGSOb<ZT,FT>::u;
  using MatGSOb<ZT,FT>::enable_transform;
  using MatGSOb<ZT,FT>::cols_locked; // maybe scratch.
  using MatGSOb<ZT,FT>::enable_int_gram;
  using MatGSOb<ZT,FT>::gso_valid_cols;
  using MatGSOb<ZT,FT>::enable_inverse_transform;
  using MatGSOb<ZT,FT>::u_inv_t;
  using MatGSOb<ZT,FT>::sym_g;
  using MatGSOb<ZT,FT>::mu;
  using MatGSOb<ZT,FT>::r;
  using MatGSOb<ZT,FT>::ztmp1;
  using MatGSOb<ZT,FT>::ztmp2;
  using MatGSOb<ZT,FT>::row_op_force_long;
  using MatGSOb<ZT,FT>::alloc_dim;
  using MatGSOb<ZT,FT>::get_mu;
  using MatGSOb<ZT,FT>::get_r;
  using MatGSOb<ZT,FT>::gptr;  
  using MatGSOb<ZT,FT>::invalidate_gso_row;

  using MatGSOb<ZT,FT>::create_row;
  using MatGSOb<ZT,FT>::remove_last_row;
  using MatGSOb<ZT,FT>::print_mu_r_g;
  using MatGSOb<ZT,FT>::discover_all_rows;
  using MatGSOb<ZT,FT>::update_gso;
  using MatGSOb<ZT,FT>::update_gso_row;
  using MatGSOb<ZT,FT>::symmetrize_g;


  MatGGSO(Matrix<ZT> &arg_g, Matrix<ZT> &arg_u, Matrix<ZT> &arg_uinv_t, int flags)
      : MatGSOb<ZT,FT>(arg_u,arg_uinv_t, GSO_INT_GRAM )    
  {
    FPLLL_DEBUG_CHECK(!(enable_int_gram && enable_row_expo));

    gptr = &arg_g;
    if (gptr == nullptr){ cerr << "Pointer of the Grammatrix equal to the nullpointer.\n"; exit(1); } 
    d = gptr->get_rows();
    size_increased();

#ifdef DEBUG
    row_op_first = row_op_last = -1;
#endif
  }


public:



  virtual inline long  get_max_exp_of_b();
  virtual inline bool b_row_is_zero(int i);
  virtual inline int get_cols_of_b();
  virtual inline int get_rows_of_b();
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

  virtual inline FT &get_gram(FT &f, int i, int j);
  
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
  // b[i] <-> b[j] (i < j)
  virtual void row_swap(int i, int j);
  //virtual void apply_transform(const Matrix<FT> &transform, int src_base, int target_base);

};
  

template <class ZT, class FT> inline long MatGGSO<ZT,FT>::get_max_exp_of_b()
{
  if (gptr == nullptr) { cerr << "Error: gptr is equal to the nullpointer.\n"; exit(1); }
  Matrix<ZT> &g = *gptr;
  // normally, this returns the maximum
  // exponent of b. We approximate this temporarily 
  // by the half of the maximum exponent of g.
  // beware of errors, because maybe g is
  // zero in the upper half of the matrix.
  return g.get_max_exp() / 2;
}
template <class ZT, class FT> inline bool MatGGSO<ZT,FT>::b_row_is_zero(int i)
{
  if (gptr == nullptr) { cerr << "Error: gptr is equal to the nullpointer.\n"; exit(1); }
  Matrix<ZT> &g = *gptr;
  // normally this returns whether the
  // ith row of the basis is the zero row.
  // now, we just check whether g[i][i] is zero.
  // beware of errors, because maybe
  // the function is_zero is not cdefined
  // for type of g[i][i]
  return g[i][i].is_zero();
}
template <class ZT, class FT> inline int MatGGSO<ZT,FT>::get_cols_of_b()
{
  if (gptr == nullptr) { cerr << "Error: gptr is equal to the nullpointer.\n"; exit(1); }
  // in this gram-matrix version, the number
  // of columns of b is the same as the
  // number of colums of g.
  return gptr->get_cols();
}

template <class ZT, class FT> inline int MatGGSO<ZT,FT>::get_rows_of_b()
{
  if (gptr == nullptr) { cerr << "Error: gptr is equal to the nullpointer.\n"; exit(1); }
  // in this gram-matrix version, the number
  // of columns of b is the same as the
  // number of colums of g.
  return gptr->get_rows();
}

template <class ZT, class FT> inline void MatGGSO<ZT,FT>::negate_row_of_b(int i)
{
      if (enable_int_gram) {
        
        for (int j = 0; j < get_rows_of_b(); j++) {
          if (j != i) {
            sym_g(i,j).neg(sym_g(i,j));
          }
        } 
      }

}

template <class ZT, class FT> inline FT &MatGGSO<ZT, FT>::get_gram(FT &f, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && j >= 0 && j <= i && j < n_source_rows &&
                    !in_row_op_range(i));
  if (enable_int_gram) 
  {
    if (gptr == nullptr) {
        cerr << "Error: gptr is equal to the nullpointer."; exit(1); 
    }
    f.set_z((*gptr)(i, j));
  }
  return f;
}

template <class ZT, class FT> inline void MatGGSO<ZT, FT>::remove_last_rows(int n_removed_rows)
{
  FPLLL_DEBUG_CHECK(!cols_locked && d >= n_removed_rows);
  d -= n_removed_rows;
  n_known_rows  = min(n_known_rows, d);
  n_source_rows = n_known_rows;
  if (enable_transform)
    u.set_rows(d);
}

template <class ZT, class FT> inline void MatGGSO<ZT, FT>::create_rows(int n_new_rows)
{
  FPLLL_DEBUG_CHECK(!cols_locked);
  int old_d = d;
  d += n_new_rows;

  if (enable_transform)
  {
    u.set_rows(d);
    for (int i = old_d; i < d; i++)
      for (int j = 0; j < u.get_cols(); j++)
        u[i][j]  = 0;
  }
  size_increased(); 
  if (n_known_rows == old_d)
    MatGSOb<ZT,FT>::discover_all_rows();
}

FPLLL_END_NAMESPACE

#endif
