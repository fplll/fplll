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

#include "gso_givens.h"

FPLLL_BEGIN_NAMESPACE

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::add_operation(FT c, FT s, int row, int col)
{
  ops_c(row, col) = c;
  ops_s(row, col) = s;
}

// The givens rotation introduces a zero in place  (k,j) in the matrix, and puts the residue
// on place (k,i). Most of the cases, j = i-1
//              j|i
//   _________________________
//  |                         |
//  |_________________________|
// k |__________|x|y|__________|
//  |                         |
//  |                         |
//  |_________________________|
//
//        TRAMSFORMS TO
//
//
//              j|i
//   _________________________
//  |           * *           |
//  |_________________________|
// k |__________|x|0|__________|
//  |           * *           |
//  |           * *           |
//  |_________________________|
//
//
//
//  The stars indicate that the values
//  there can be changed [in linear combinations].
// The rest of  the matrix is fixed.
//

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::givens_rotation(int row, int col)
{
  if (l_givens(row, col).is_zero())
  {
    add_operation(1.0, 0.0, row, col);
    return;
  }

  FT c, s, tempvar;

  // computes c,s
  tempvar.hypot(l_givens(row, col - 1), l_givens(row, col));
  c.div(l_givens(row, col - 1), tempvar);
  s.div(l_givens(row, col), tempvar);

  // applies the c,s-rotation to place (row,col)
  virtual_givens_rotation(row, col, c, s);

  // adds the operation to the operation-matrices.
  add_operation(c, s, row, col);

  // force zero on the place where a zero should be.
  l_givens(row, col) = 0.0;
}

// Applies a 'virtual' givens rotation, meaning
// that the c and s are not from the row to which
// you apply the rotation.
template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::virtual_givens_rotation(int row, int col, FT c, FT s)
{
  FT tempvar, tempvar2;
  tempvar  = l_givens(row, col - 1);
  tempvar2 = l_givens(row, col);
  l_givens(row, col - 1).mul(tempvar, c);
  l_givens(row, col - 1).addmul(tempvar2, s);

  l_givens(row, col).neg(s);
  l_givens(row, col).mul(tempvar, l_givens(row, col));
  l_givens(row, col).addmul(tempvar2, c);
}

// The givens_row_reduction introduces a zero-sequence from the diagonal of row_k
// to the 'rightmost_nonzero_entry' column.
//
//              k|        r
//   _________________________
//  |                         |
//  |_________________________|
// k |__________|a|b|c|d|e|f|__|
//  |                         |
//  |                         |
//  |_________________________|
//
//        TRAMSFORMS TO
//
//
//              k|        r
//   _________________________
//  |           * * * * * *   |
//  |_________________________|
// k |__________|z|0|0|0|0|0|__|
//  |           * * * * * *   |
//  |           * * * * * *   |
//  |-------------------------|
//
// The stars indicate that the entries
// there might be changed into linear
// combinations of their neighbour-entries.
// The rest of the matrix is fixed.

template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::givens_row(int row, int rightmost_nonzero_entry)
{
  for (int i = rightmost_nonzero_entry; i > row; i--)
    givens_rotation(row, i);
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::givens_row(int row)
{
  givens_row(row, l_givens.get_cols() - 1);
}

// Applies to row 'row' all givens operation
// from the previous rows. This is because
// we postpone the column operations up to the
// moment that we really need them.

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::apply_givens_operations(int row)
{

  for (int i = 0; i < row; i++)
  {
    for (int j = l_givens.get_cols() - 1; j > i; j--)
    {
      // maybe do a check if ops_s(i,j) is zero?
      virtual_givens_rotation(row, j, ops_c(i, j), ops_s(i, j));
    }
  }
}

/*
template <class ZT, class FT> void MatGSOGivens<ZT,FT>::apply_givens_operations(int start_row, int
end_row){
        FT tempvar,tempvar2;
    for(int i = 0; i < start_row; i++) {
    for(int j = l_givens.get_cols()-1; j>i; j--) {
      // maybe do a check if ops_s(i,j) is zero?

      for(int k = start_row; k < end_row; k++) {
        tempvar = l_givens(k, j-1);
        tempvar2 = l_givens(k, j);
        l_givens(k, j-1).mul(tempvar, ops_c(i,j));
        l_givens(k, j-1).addmul(tempvar2, ops_s(i,j));

        l_givens(k, j).neg(ops_s(i,j));
        l_givens(k, j).mul(tempvar, l_givens(k, j));
        l_givens(k, j).addmul(tempvar2, ops_c(i,j));
      }
    }
    }


}
*/
/*
template <class ZT, class FT> void MatGSOGivens<ZT,FT>::triangularize(int start_row, int end_row){
        FT tempvar,tempvar2,c,s;

    for(int i = start_row; i < end_row; i++) {
        for(int j = l_givens.get_cols()-1; j>i; j--) {
          // maybe do a check if ops_s(i,j) is zero?
          if (l_givens(i,j).is_zero()) {
            add_operation(1.0,0.0,i,j);
            continue;
          }
          tempvar.hypot(l_givens(i, j-1), l_givens(i, j));
          ops_c(i,j).div(l_givens(i, j-1), tempvar);
          ops_s(i,j).div(l_givens(i, j), tempvar);
          c = ops_c(i,j);
          s = ops_s(i,j);

          for(int k = i; k < end_row; k++) {
            tempvar = l_givens(k, j-1);
            tempvar2 = l_givens(k, j);
            l_givens(k, j-1).mul(tempvar, c);
            l_givens(k, j-1).addmul(tempvar2, s);

            l_givens(k, j).neg(s);
            l_givens(k, j).mul(tempvar, l_givens(k, j));
            l_givens(k, j).addmul(tempvar2, c);
          }
          l_givens(i,j) = 0.0; // force 0.
        }
    compute_mu_and_r(i);
    }

}
*/

template <class ZT, class FT> bool MatGSOGivens<ZT, FT>::real_update_gso_row(int row)
{
  // copies b[row] to l_givens[row]
  copy_b_to_l_givens(row);

  // applies all previous givens operations to l_givens[row]
  apply_givens_operations(row);

  // givens-reduces l_givens[row], and stores (c,s) values in the ops-matrices
  givens_row(row);

  // computes mu and r.
  // TODO to make lazy!!!
  compute_mu_and_r(row);

  // Discover the new row.
  if (row >= n_known_rows)
  {
    discover_row();
  }
  return true;
}

template <class ZT, class FT> bool MatGSOGivens<ZT, FT>::update_gso_row(int row, int last_j)
{

  if (full_lazy)
  {
    if (n_known_rows < d)
    {  // Only the first time the full GSO is computed
      for (int i = n_known_rows; i < d; i++)
        real_update_gso_row(i);
    }
    return true;  // Rest of the time, do NO recomputation.
  }

  if (move_lazy)
  {
    // If you did 'move´ lazy, the Givens-operations of earlier rows need to be recomputed...!
    for (int i = lazy_row_start; i < row; i++)
      real_update_gso_row(i);
    lazy_row_start = d;
  }

  return real_update_gso_row(row);
}

template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::recompute_givens_matrix(int start_row, int last_row)
{
  for (int i = start_row; i < last_row; i++)
  {
    real_update_gso_row(i);
  }
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::recompute_givens_matrix(int last_row)
{
  recompute_givens_matrix(0, last_row);
}

// Includes also row-exponent compatibility.
// Essentially copies b[row] to l_givens[row].

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::copy_b_to_l_givens(int row)
{
  if (enable_row_expo)
  {
    int n         = b.get_cols();
    long max_expo = LONG_MIN;
    for (int j = 0; j < n; j++)
    {
      b(row, j).get_f_exp(l_givens(row, j), tmp_col_expo[j]);
      max_expo = max(max_expo, tmp_col_expo[j]);
    }
    for (int j = 0; j < n; j++)
    {
      l_givens(row, j).mul_2si(l_givens(row, j), tmp_col_expo[j] - max_expo);
    }
    row_expo[row] = max_expo;

    // normalize_givens_row(row);
  }
  else
  {
    for (int j = 0; j < l_givens.get_cols(); j++)
    {
      l_givens(row, j).set_z(b(row, j));
    }
  }
}

template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::compute_mu_and_r_columns(int starting_column, int last_column)
{
  for (int col = starting_column; col <= last_column; col++)
  {
    r(col, col).mul(l_givens(col, col), l_givens(col, col));
    for (int i = col; i < n_known_rows; i++)
      mu(i, col).div(l_givens(i, col), l_givens(col, col));

    ftmp1 = l_givens(col, col);
    for (int i = col; i < n_known_rows; i++)
      r(i, col).mul(ftmp1, l_givens(i, col));
  }
}

// TODO MAKE LAZY

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::compute_mu_and_r(int row)
{
  for (int k = 0; k < row; k++)
  {
    mu[row][k].div(l_givens[row][k], l_givens(k, k));
    r[row][k].mul(l_givens[row][k], l_givens(k, k));
  }
  r[row][row].mul(l_givens[row][row], l_givens[row][row]);
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::update_bf(int i)
{

  if (enable_row_expo)
  {
    int n = b.get_cols();

    // long max_expo = LONG_MIN;
    for (int j = 0; j < n; j++)
    {
      b(i, j).get_f_exp(bf(i, j), tmp_col_expo[j]);
      // max_expo = max(max_expo, tmp_col_expo[j]);
    }
    for (int j = 0; j < b.get_cols(); j++)
      bf(i, j).mul_2si(bf(i, j), tmp_col_expo[j] - row_expo[i]);
    //{
    //  bf(i, j).mul_2si(bf(i, j), - );
    //}
    // row_expo[i] = max_expo;
  }
  else
  {
    for (int j = 0; j < b.get_cols(); j++)
    {
      bf(i, j).set_z(b(i, j));
    }
  }
}

// "Full columns" ... does it improve time?

template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::full_column_givens_rotation(int row, int col)
{

  FT c, s, tempvar;
  tempvar.hypot(l_givens(row, col - 1), l_givens(row, col));
  c.div(l_givens(row, col - 1), tempvar);
  s.div(l_givens(row, col), tempvar);

  for (int k = row; k < n_known_rows; k++)
  {
    ftmp1 = l_givens(k, col - 1);
    ftmp2 = l_givens(k, col);
    l_givens(k, col - 1).mul(ftmp1, c);  // r_(k,col_i) = c*r_(k,col_i) + s*r_(k,col_j)
    l_givens(k, col - 1).addmul(ftmp2, s);

    l_givens(k, col).neg(s);
    l_givens(k, col).mul(ftmp1, l_givens(k, col));
    l_givens(k, col).addmul(ftmp2, c);  // r_(k,col_j) = -s*r_(k,col_i) + c*r_(k,col_j)
  }
  // "Forcing" zero
  // l_givens(row, col) = 0.0;
}

template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::full_column_givens_row(int row, int rightmost_nonzero_entry)
{
  for (int i = rightmost_nonzero_entry; i > row; i--)
    full_column_givens_rotation(row, i);
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::invalidate_gram_row(int i)
{
  // TODO maybe some functionality later.
}

// Givens ready.
template <class ZT, class FT> void MatGSOGivens<ZT, FT>::row_add(int i, int j)
{

  b[i].add(b[j], b.get_cols());
  if (enable_transform)
  {
    u[i].add(u[j]);
    if (enable_inverse_transform)
      u_inv_t[j].sub(u_inv_t[i]);
  }

  if (j < i)
  {
    if (enable_row_expo)
    {

      l_givens[i].addmul_2exp(l_givens[j], 1.0, row_expo[j] - row_expo[i], ftmp1);

      // TODO LAZY
      compute_mu_and_r(i);
    }
    else
    {
      // Doing b_i <- b_i + c b_j doesn't affect the triangularity

      l_givens[i].add(l_givens[j], j + 1);

      // TODO: making this lazy
      compute_mu_and_r(i);
    }
  }
  else
  {
    throw std::runtime_error("Error: i < j in row_add");
  }
}

// Givens ready
template <class ZT, class FT> void MatGSOGivens<ZT, FT>::row_sub(int i, int j)
{

  b[i].sub(b[j], b.get_cols());
  if (enable_transform)
  {
    u[i].sub(u[j]);
    if (enable_inverse_transform)
      u_inv_t[j].add(u_inv_t[i]);
  }

  if (j < i)
  {
    if (enable_row_expo)
    {
      l_givens[i].addmul_2exp(l_givens[j], -1.0, row_expo[j] - row_expo[i], ftmp1);
      compute_mu_and_r(i);
    }
    else
    {
      // Doing b_i <- b_i + c b_j doesn't affect the triangularity
      l_givens[i].sub(l_givens[j], j + 1);

      // TODO: making this lazy

      compute_mu_and_r(i);
    }
  }
  else
  {
    // i < j, so affects triangularity
    throw std::runtime_error("Error: i < j in row_sub");
  }
}

// Givens-ready
template <class ZT, class FT> void MatGSOGivens<ZT, FT>::row_addmul_si(int i, int j, long x)
{

  // TODO addmul_si not possible, because
  // l_givens is of type NumVect
  // and not of MatrixRow
  // Is this well-resolved??

  b[i].addmul_si(b[j], x, b.get_cols());

  if (enable_transform)
  {
    u[i].addmul_si(u[j], x);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si(u_inv_t[i], -x);
  }

  if (j < i)
  {
    if (enable_row_expo)
    {
      l_givens[i].addmul_2exp(l_givens[j], x, row_expo[j] - row_expo[i], ftmp1);
      compute_mu_and_r(i);
    }
    else
    {

      // Doing b_i <- b_i + c b_j doesn't affect the

      // TODO Changing to addmul okay?
      l_givens[i].addmul(l_givens[j], x, j + 1);  // j+1 really needed!

      // TODO: making this lazy
      compute_mu_and_r(i);
    }
  }
  else
  {
    throw std::runtime_error("Error: i < j in row_addmul_si");
  }
}

// TODO Still needs to be implemented!
template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::row_addmul_si_2exp(int i, int j, long x, long expo)
{
  throw std::runtime_error(
      "Error in row_addmul_si_2exp: exponents are not yet implemented for givens rotations");
  b[i].addmul_si_2exp(b[j], x, expo, b.get_cols(), ztmp1);
  if (enable_transform)
  {
    u[i].addmul_si_2exp(u[j], x, expo, ztmp1);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si_2exp(u_inv_t[i], -x, expo, ztmp1);
  }

  if (j < i)
  {
    if (enable_row_expo)
    {
      long double dx = (long double)x;

      l_givens[i].addmul_2exp(l_givens[j], dx, row_expo[j] - row_expo[i] + expo, ftmp1);
      compute_mu_and_r(i);
    }
    else
    {
      // Doing b_i <- b_i + c b_j doesn't affect the triangularity

      l_givens[i].addmul_2exp(l_givens[j], x, expo, ftmp1);

      // TODO: making this lazy
      compute_mu_and_r(i);
    }
  }
  else
  {
    throw std::runtime_error("Error: i < j in row_addmul_si");
  }
}

template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::row_addmul_2exp(int i, int j, const ZT &x, long expo)
{
  b[i].addmul_2exp(b[j], x, expo, b.get_cols(), ztmp1);
  if (enable_transform)
  {
    u[i].addmul_2exp(u[j], x, expo, ztmp1);
    if (enable_inverse_transform)
    {
      ZT minus_x;
      minus_x.neg(x);
      u_inv_t[j].addmul_2exp(u_inv_t[i], minus_x, expo, ztmp1);
    }
  }

  FT tmpx;
  tmpx.set_z(x);

  if (j < i)
  {
    // Doing b_i <- b_i + c b_j doesn't affect the triangularity

    if (enable_row_expo)
    {
      long double dx = x.get_d();

      l_givens[i].addmul_2exp(l_givens[j], dx, row_expo[j] - row_expo[i] + expo, ftmp1);
      // LAZY
      compute_mu_and_r(i);
    }
    else
    {

      for (int k = 0; k < j; k++)
      {
        ftmp1.mul(l_givens(j, k), tmpx);
        ftmp1.mul_2si(ftmp1, expo);
        l_givens(i, k).add(l_givens(i, k), ftmp1);
      }
      // LAZY
      compute_mu_and_r(i);
    }
  }
  else
  {
    throw std::runtime_error("Error: i < j in row_addmul_2exp");
  }
}

// in row_addmul_we we must have i > j
template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::row_addmul_we(int i, int j, const FT &x, long expo_add)
{

  FPLLL_DEBUG_CHECK(j >= 0 && /* i > j &&*/ i < n_known_rows && j < n_source_rows);
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

// In row_swap, i < j
// Is Givens-ready.
template <class ZT, class FT> void MatGSOGivens<ZT, FT>::row_swap(int i, int j)
{
  FPLLL_DEBUG_CHECK(!enable_inverse_transform);
  if (j < i)
  {
    int k = i;
    i     = j;
    j     = k;
  }

  // *******************
  // Swap the rows of b
  // ********************
  b.swap_rows(i, j);
  if (enable_transform)
  {
    u.swap_rows(i, j);
  }
  // *****************
  // Givens equivalent
  // *****************
  l_givens.swap_rows(i, j);
  mu.swap_rows(i, j);
  r.swap_rows(i, j);

  if (move_lazy)
  {

    full_column_givens_row(i, j);
    for (int k = i + 1; k < j; k++)
    {
      full_column_givens_rotation(k, k + 1);
    }
    compute_mu_and_r_columns(i, j);
  }
  else
  {
    recompute_givens_matrix(i, j);
  }

  // if (always_recompute){
  //  recompute_givens_matrix();
  //}
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::move_row(int old_r, int new_r)
{
  // cerr << "Move " << old_r << " to " << new_r << " and n_known_rows = " << n_known_rows << endl;
  FPLLL_DEBUG_CHECK(!cols_locked);
  if (new_r < old_r)
  {
    FPLLL_DEBUG_CHECK(old_r < n_known_rows);
    b.rotate_right(new_r, old_r);

    if (enable_transform)
    {
      u.rotate_right(new_r, old_r);
      if (enable_inverse_transform)
        u_inv_t.rotate_right(new_r, old_r);
    }

    l_givens.rotate_right(new_r, old_r);

    mu.rotate_right(new_r, old_r);
    r.rotate_right(new_r, old_r);

    if (enable_row_expo)
    {
      rotate(row_expo.begin() + new_r, row_expo.begin() + old_r, row_expo.begin() + old_r + 1);
    }

    if (move_lazy)
    {
      full_column_givens_row(new_r, old_r);
      compute_mu_and_r_columns(new_r, old_r);
      lazy_row_start = min(new_r, lazy_row_start);
    }
    else
    {
      recompute_givens_matrix(new_r, n_known_rows);
    }
  }
  else if (new_r > old_r)
  {
    // throw std::runtime_error("Wrong order of rotation!");

    b.rotate_left(old_r, new_r);
    if (enable_transform)
    {
      u.rotate_left(old_r, new_r);
      if (enable_inverse_transform)
        u_inv_t.rotate_left(old_r, new_r);
    }

    l_givens.rotate_left(old_r, new_r);
    mu.rotate_left(old_r, new_r);
    r.rotate_left(old_r, new_r);
    if (enable_row_expo)
    {
      rotate(row_expo.begin() + old_r, row_expo.begin() + old_r + 1, row_expo.begin() + new_r + 1);
    }

    if (move_lazy)
    {
      for (int i = old_r; i < new_r; i++)
        full_column_givens_row(i, i + 1);

      compute_mu_and_r_columns(old_r, new_r);
      lazy_row_start = min(old_r, lazy_row_start);
    }
    else
    {
      recompute_givens_matrix(old_r, n_known_rows);
    }
  }

  if (new_r >= n_known_rows)
  {
    if (old_r < n_known_rows)
    {
      n_known_rows--;
    }
  }
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::size_increased()
{

  if (d > alloc_dim)
  {

    mu.resize(d, b.get_cols());
    l_givens.resize(d, b.get_cols());
    r.resize(d, b.get_cols());
    bf.resize(b.get_rows(), b.get_cols());
    gf.resize(b.get_rows(), b.get_rows());
    if (enable_row_expo)
    {
      row_expo.resize(d);
    }
    ops_s.resize(d, b.get_cols());
    ops_c.resize(d, b.get_cols());
    // lazy_row_c.resize(1,b.get_cols());
    // lazy_row_s.resize(1,b.get_cols());
    alloc_dim = d;
  }
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::row_op_end(int first, int last)
{
#ifdef DEBUG
  FPLLL_DEBUG_CHECK(row_op_first == first && row_op_last == last);
  row_op_first = row_op_last = -1;
#endif
  if (always_recompute)
  {
    recompute_givens_matrix(first, last);
  }
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::discover_row()
{
  FPLLL_DEBUG_CHECK(n_known_rows < d);
  FPLLL_DEBUG_CHECK(!(cols_locked));
  n_known_rows++;
}

template class MatGSOGivens<Z_NR<long>, FP_NR<double>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<double>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class MatGSOGivens<Z_NR<long>, FP_NR<long double>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<long double>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<long double>>;

#endif

#ifdef FPLLL_WITH_QD
template class MatGSOGivens<Z_NR<long>, FP_NR<dd_real>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<dd_real>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<dd_real>>;

template class MatGSOGivens<Z_NR<long>, FP_NR<qd_real>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<qd_real>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class MatGSOGivens<Z_NR<long>, FP_NR<dpe_t>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<dpe_t>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<dpe_t>>;
#endif

template class MatGSOGivens<Z_NR<long>, FP_NR<mpfr_t>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<mpfr_t>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<mpfr_t>>;

FPLLL_END_NAMESPACE
