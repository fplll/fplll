/*
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

/* source file */
/* to be templated one day if the need arises */

#include "hnf.h"
#include "util.h"

FPLLL_BEGIN_NAMESPACE

// the following macro is a loop to navigate inside a matrix in a diagonal manner.
// INDEXES :
//  - j is a column index, l the minimum value of j, k the pivot row index.
//  - r,c must be the number of rows and column of the input matrix
//  - pivot coefficient is at line k column j
//  - navigation goes diagonally backwards from the last row last column
// INSTRUCTIONS :
//  - PRE_INST : instructions before checking the pivot value
//  - NEG_INST : instructions if pivot value is negative
//  - ZERO_INST : instructions if pivot value is zero
//  - POS_INST : instructions if pivot value is positive

#define _DIAGONAL_INDEX_LOOP_(B, PRE_INST, NEG_INST, ZERO_INST, POS_INST)                          \
  {                                                                                                \
    for (j = c - 1, k = r - 1; j >= l; j--, k--)                                                   \
    {                                                                                              \
      /* operations on pivot position before any value check */                                    \
      {                                                                                            \
        PRE_INST                                                                                   \
      }                                                                                            \
      /* checks if diagonal is non-negative */                                                     \
      if (B[k][j] < 0)                                                                             \
      {                                                                                            \
        {                                                                                          \
          NEG_INST                                                                                 \
        }                                                                                          \
      }                                                                                            \
      /* if "diagonal" entry is zero */                                                            \
      else if (B[k][j] == 0)                                                                       \
      {                                                                                            \
        {                                                                                          \
          ZERO_INST                                                                                \
        }                                                                                          \
      }                                                                                            \
      else                                                                                         \
      /* if "diagonal" entry is positive */                                                        \
      {                                                                                            \
        {                                                                                          \
          POS_INST                                                                                 \
        }                                                                                          \
      }                                                                                            \
    }                                                                                              \
  }

/* clang-format off */
// HNF checking macro to check if B is a HNF, and add an instruction if the pivot is correct
#define _HNF_CHECK_(B, INSTRUCTION_AT_PIVOT)                                                       \
  {                                                                                                \
    _DIAGONAL_INDEX_LOOP_(B                                                                        \
      ,                                                                                            \
      /* checks if everything above diagonal is zero */                                            \
      for (i = k - 1; i >= 0; i--) {                                                               \
        if (B[i][j] != 0) { status = RED_HNF_FAILURE;}                                                             \
      }                                                                                            \
      ,                                                                                            \
      /* neg : diagonal shouldn't be negative */                                                   \
      { status = RED_HNF_FAILURE; }                                                                                \
      ,                                                                                            \
      /* zero : pivot row position doesn't change, increment to compensate */                      \
      /* lower the limit as well, we skipped one column without increasing row */                  \
      { k++; if (l > 0) { l--; } }                                                                 \
      ,                                                                                            \
      /* pos : checks if everything below is also positive and lower than the pivot */             \
      /* leave a potential extra instruction */                                                    \
      for (i = k + 1; i < r; i++) {                                                                \
        if ((B[i][j] > B[k][j]) || (B[i][j] < 0)) { status = RED_HNF_FAILURE; }                                    \
      }                                                                                            \
      /* leave a potential extra instruction */                                                    \
      {INSTRUCTION_AT_PIVOT}                                                                       \
    )                                                                                              \
  }

/* clang-format on */

// HNF computation macro, reduce coefficients of row by the pivot starting from column j
#define _REDUCE_ROW_(row, j, j2, pivot)                                                            \
  {                                                                                                \
    q.fdiv_q(row[j], pivot[j]);                                                                    \
    for (j2 = j; j2 >= 0; j2--)                                                                    \
    {                                                                                              \
      row[j2].submul(q, pivot[j2]);                                                                \
    }                                                                                              \
  }

// HNF computation macro, reduce lower coefficients of column j by the pivot in row k
#define _REDUCE_LOWER_(B, i, j, j2, k, r)                                                          \
  {                                                                                                \
    for (i = k + 1; i < r; i++)                                                                    \
    {                                                                                              \
      _REDUCE_ROW_(B[i], j, j2, B[k])                                                              \
    }                                                                                              \
  }

// Takes a matrix and two row indexes and applies the xgcd reduction, and updates pivot
// j is the head coefficient's column. d,u,v,r1d,r2d are variables used for xgcd
#define _REDUCE_VECT_XGCD(B, to_red, j, j2, pivot, b, d, u, v, r1d, r2d)                           \
  {                                                                                                \
    /* skip zero */                                                                                \
    if (B[to_red][j] != 0)                                                                         \
    {                                                                                              \
      /* reduce row to_red with row pivot */                                                       \
      d.xgcd(u, v, B[pivot][j], B[to_red][j]);                                                     \
      /* This is in case the pivot is the gcd, saves some time (useful for ones for sure) */       \
      if (d.cmpabs(B[pivot][j]) == 0)                                                              \
      {                                                                                            \
        b.divexact(B[to_red][j], B[pivot][j]);                                                     \
        for (j2 = j; j2 >= 0; j2--)                                                                \
        {                                                                                          \
          B[to_red][j2].submul(b, B[pivot][j2]);                                                   \
        }                                                                                          \
      }                                                                                            \
      else                                                                                         \
      {                                                                                            \
        r2d.divexact(B[to_red][j], d);                                                             \
        r1d.divexact(B[pivot][j], d);                                                              \
        /* only compute relevant values (rest should be guaranteed zero) */                        \
        for (j2 = j; j2 >= 0; j2--)                                                                \
        {                                                                                          \
          /* precompute values */                                                                  \
          b.mul(u, B[pivot][j2]);                                                                  \
          b.addmul(v, B[to_red][j2]);                                                              \
          /* new vector i-1 value */                                                               \
          B[to_red][j2].mul(r1d, B[to_red][j2]);                                                   \
          B[to_red][j2].submul(r2d, B[pivot][j2]);                                                 \
          /* new vector i value */                                                                 \
          B[pivot][j2] = b;                                                                        \
        }                                                                                          \
      }                                                                                            \
    }                                                                                              \
  }

// Class method templates
template <class ZT> bool HNFReduction<ZT>::is_reduced()
{
  status = RED_SUCCESS;
  /* matrix bounds */
  int r = basis.get_rows(), c = basis.get_cols();
  /* matrix indexes (k = pivot)*/
  int i, j, k;
  /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);
  _HNF_CHECK_(basis, {});
  return status == RED_SUCCESS;
}

template <class ZT> void HNFReduction<ZT>::hnf()
{
  status = RED_HNF_FAILURE;
  if (method == HM_AUTO)
  {
    status = hnf_autoselect(basis);
  }
  else if (method == HM_XGCD)
  {
    status = hnf_xgcd_reduction(basis);
  }
  else if (method == HM_CLASSIC)
  {
    status = hnf_classical_reduction(basis);
  }
  else if (method == HM_MINORS)
  {
    status = hnf_minors_reduction(basis);
  }
  else if (method == HM_MODULO)
  {
    status = hnf_modular_reduction(basis, det);
  }
  // rest in development
  else if (method == HM_PERNETSTEIN)
  {
    cerr << "HNF method not implemented yet\n";
  }
  // not existant
  else
  {
    cerr << "invalid method\n";
  }
}

template <class ZT> bool HNFReduction<ZT>::in_lattice_given_hnf(const vector<Z_NR<ZT>> &w)
{
  status = RED_SUCCESS;
  /* matrix bounds */
  int r = basis.get_rows(), c = basis.get_cols();
  /* matrix indexes (k = pivot)*/
  int i, j, j2, k;
  /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);

  vector<Z_NR<ZT>> v = w;
  Z_NR<ZT> q;

  /* test if the vector size is correct */
  if ((int)v.size() != r)
  {
    cerr << "in_hnf error : matrix-vector sizes do not match\n";
    status = RED_HNF_FAILURE;
  }
  /* clang-format off */
  _HNF_CHECK_(basis,
    /* reduces the vector by the row pivot vector for membership test */
    _REDUCE_ROW_(v, j, j2, basis[k])
  );
  /* clang-format on */

  // membership test : the vector after reduction should be a zero one.
  for (j = 0; j < c; j++)
  {
    if (v[j] != 0)
    {
      status = RED_HNF_FAILURE;
    }
  }

  return status == RED_SUCCESS;
}

template <class ZT> bool HNFReduction<ZT>::in_lattice_given_hnf(const ZZ_mat<ZT> &A)
{
  status = RED_SUCCESS;

  ZZ_mat<ZT> tmp = A;
  Z_NR<ZT> q;

  /* matrix bounds */
  int r = basis.get_rows(), c = basis.get_cols();
  /* matrix indexes (k = pivot)*/
  int i, j, j2, k;
  /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);

  /* test if the matrix size is correct */
  if ((A.get_cols() != basis.get_cols()) || (A.get_rows() != basis.get_rows()))
  {
    cerr << "in_hnf error : matrix-matrix sizes do not match\n";
    status = RED_HNF_FAILURE;
  }

  /* clang-format off */
  _HNF_CHECK_(basis,
    /* reduces the matrix tmp by the row pivot vector of A for membership test */
    for (i = 0; i < r; i++)
      {
       _REDUCE_ROW_(tmp[i], j, j2, basis[k])
      }
      // cout << "HNF reduction for pivot " << k << endl << B << endl;
      // cout << "tmp reduction for pivot " << k << endl << tmp << endl;
  );
  /* clang-format on */

  // membership test : the matrix after reduction should be a zero one.
  for (i = 0; i < r; i++)
  {
    for (j = 0; j < c; j++)
    {
      if (tmp[i][j] != 0)
      {
        status = RED_HNF_FAILURE;
      }
    }
  }

  return status == RED_SUCCESS;
}

template <class ZT> int hnf_xgcd_reduction(ZZ_mat<ZT> &B)
{
  /* matrix indexes (k = pivot)*/
  int i, j, j2, k;
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols();
  /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);
  /* ZT integers for operations */
  Z_NR<ZT> r1d, r2d, b, u, v, d, q;

  /* initializes the transformation matrix */
  // ZZ_mat<ZT> U;
  // U = B;
  // U.gen_identity(B.get_rows());

  /* clang-format off */
  _DIAGONAL_INDEX_LOOP_(B,
    /* pre-check : iterates above the pivot to construct it */
    for (i = k - 1; i >= 0; i--) {
      /* reduce row i + 1 with row i and makes i the non zero vector*/
      _REDUCE_VECT_XGCD(B,(i+1),j,j2,(i),b,d,u,v,r1d,r2d)
    }
    /* swap first row with the pivot row */
    B.swap_rows(0, k);
    ,
    /* neg : change sign of the row vector if the diagonal entry is negative */
    for (j2 = j; j2 >= 0; j2--) { B[k][j2].neg(B[k][j2]); }
    _REDUCE_LOWER_(B, i, j, j2, k, r)
    ,
    /* zero : modify pivot position */
    { k++; if (l > 0) { l--; } }
    ,
    /* pos : reduce lower entries of column j with row k */
    _REDUCE_LOWER_(B, i, j, j2, k, r)
  )
  /* clang-format on */

  // return in_lattice_given_hnf(B,U);
  return 0;
}

template <class ZT> int hnf_classical_reduction(ZZ_mat<ZT> &B)
{
  /* matrix indexes (k = pivot)*/
  int i, j, j2, k, i_min;
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols();
  /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);

  /* ZT integers for operations */
  Z_NR<ZT> q;
  /* value to check whether or not a column is already reduced */
  int reduced;

  /* initializes the transformation matrix */
  // ZZ_mat<ZT> U;
  // U = B;
  // U.gen_identity(B.get_rows());

  /* clang-format off */
  _DIAGONAL_INDEX_LOOP_(B,
    /* pre-check : iterates above the pivot to construct it */
    /* check if the column j is already reduced above*/
    for (i = k , reduced = 1; (i > 0) && reduced; )
      {
        i--;
        reduced = (B[i][j] == 0);
      }
    if (!reduced) {
      /* the first non-zero leading coefficient is either at k or i */
      if ( B[k][j] != 0 && B[i][j].cmpabs(B[k][j]) > 0)
        {
          i_min = k;
        } else {
          i_min = i;
        }
      /* have to find the row with the minimum non-zero absolute value */
      for (i = i - 1 ; i >= 0; i--) {
        if (B[i][j] != 0 && B[i_min][j].cmpabs(B[i][j]) > 0)
          { i_min = i; }
      }
      /* use the minimum row as pivot, reduce the rest above */
      B.swap_rows(i_min, k);
      for (i = k - 1, reduced = 1; (i >= 0); i--)
      {
        _REDUCE_ROW_(B[i], j, j2, B[k])
      }
      /* have to check again if we're finished before going on with the checks*/
      j++; k++; continue;
    }
    ,
    /* neg : change sign of the row vector if the diagonal entry is negative, then reduce */
    for (j2 = j; j2 >= 0; j2--) {
      B[k][j2].neg(B[k][j2]);
    }
    _REDUCE_LOWER_(B, i, j, j2, k, r)
    ,
    /* zero : modify pivot position */
    { k++; if (l > 0) { l--; } }
    ,
    /* pos : reduce lower entries of column j with row k */
    _REDUCE_LOWER_(B, i, j, j2, k, r)
  )
  /* clang-format on */

  // return in_lattice_given_hnf(B,U);
  return 0;
}

template <class ZT> int hnf_modular_reduction(ZZ_mat<ZT> &B, const Z_NR<ZT> D)
{
  /* matrix indexes (k = pivot)*/
  int i, j, k;
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols();
  if (r != c)
  {
    cerr << "modular method requires a square invertible matrix" << '\n';
    return -1;
  }

  /* ZT integers for operations */
  Z_NR<ZT> R = D, R2, d, u, v, r1d, r2d, b, q;

  /* initializes the transformation matrix */
  // ZZ_mat<ZT> U;
  // U = B;
  // U.gen_identity(B.get_rows());

  for (k = c - 1; k >= 0; k--)
  {
    /* Set a bound to not cross (about half the determinant, to minimize entry sizes) */
    R2.fdiv_q_2exp(R, 1);

    /* if the result is zero on the diagonal, the det remainder value is reached */
    if (B[k][k] == 0)
    {
      B[k][k] = R;
    }
    /* iterates above the pivot to construct it */
    for (i = k - 1; i >= 0; i--)
    {
      /* skip zeroes (if any) */
      if (B[i][k] == 0)
      {
        continue;
      }
      /* reduce row i with pivot k mod R */
      d.xgcd(u, v, B[k][k], B[i][k]);
      r1d.divexact(B[k][k], d);
      r2d.divexact(B[i][k], d);
      /* only compute relevant values (rest should be guaranteed zero) */
      for (j = k; j >= 0; j--)
      {
        /* precompute values */
        b.mul(u, B[k][j]);
        b.addmul(v, B[i][j]);
        /* new vector i value */
        B[i][j].mul(r1d, B[i][j]);
        B[i][j].submul(r2d, B[k][j]);
        B[i][j].mod(B[i][j], R);
        if (B[i][j] > R2)
        {
          B[i][j].sub(B[i][j], R);
        }
        /* new pivot vector value */
        B[k][j].mod(b, R);
        if (B[k][j] > R2)
        {
          B[k][j].sub(B[k][j], R);
        }
      }
    }

    /* refresh pivot values */
    d.xgcd(u, v, B[k][k], R);
    for (j = k; j >= 0; j--)
    {
      B[k][j].mul(u, B[k][j]);
      B[k][j].mod(B[k][j], R);
    }

    /*if the result is zero on the diagonal, it means the determinant remainder value is reached */
    if (B[k][k] == 0)
    {
      B[k][k] = R;
    }

    /* reduce higher entries of column k with pivot k */
    for (i = k + 1; i < r; i++)
    {
      _REDUCE_ROW_(B[i], k, j, B[k])
    }

    /* Sets the new determinant remainder */
    R.divexact(R, d);
  }

  // return in_lattice_given_hnf(B,U);
  return 0;
}

template <class ZT> int hnf_minors_reduction(ZZ_mat<ZT> &B)
{
  /* matrix indexes (p_i, p_j = pivot being reduced, k = block_size)*/
  int k, j, j2, l, block_size, max_minor, p_i, p_j;
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols();
  /* ZT integers for operations */
  if (r < c)
  {
    max_minor = r;
  }
  else
  {
    max_minor = c;
  }
  /* this variable indicates the row to swap with the pivot */
  int i_swap = 0;

  Z_NR<ZT> r1d, r2d, b, u, v, d, q;

  /* initializes the transformation matrix */
  // ZZ_mat<ZT> U;
  // U = B;
  // U.gen_identity(B.get_rows());

  /* put the kth principal minor in HNF (starting from bottom right corner)
   * but without computing the lower triangle for now (get diagonal coefficients)
   * k indicates the block size */
  for (block_size = 1; block_size <= max_minor; block_size++)
  {
    p_i = r - block_size;
    p_j = c - block_size;
    l   = ((c - block_size) * (c > block_size)) + 1;
    for (j = (c - 1), k = (r - 1); j >= l; j--, k--)
    {
      /* reduces p_i with row i */
      _REDUCE_VECT_XGCD(B, p_i, j, j2, k, b, d, u, v, r1d, r2d)
    }
    /* if the next last pivot is zero we swap row k for some other row (starting with the first) */
    if (B[p_i][p_j] == 0)
    {
      if (i_swap != p_i)
      {
        B.swap_rows(p_i, i_swap);
        i_swap++;
        block_size--;
        continue;
      }
      cerr << "error in minors reduction : bottom-right block is not full rank" << endl;
      return -1;
    }
    /* ensure B[p_i][p_j] is positive */
    if (B[p_i][p_j] < 0)
    {
      for (j = p_j; j >= 0; j--)
      {
        B[p_i][j].neg(B[p_i][j]);
      }
    }
    i_swap = 0;
  }

  l = ((c - r) * (c > r));
  /* finishing the lower triangle and extra rows */
  for (j = (c - 1), k = (r - 1); j >= l; j--, k--)
  {
    /* reduce extra rows */
    for (block_size = 0; block_size < r - max_minor; block_size++)
    {
      _REDUCE_VECT_XGCD(B, block_size, j, j2, k, b, d, u, v, r1d, r2d)
    }

    /* reduce below diagonal elements */
    _REDUCE_LOWER_(B, p_i, j, j2, k, r);
  }

  return 0;
}

template <class ZT> void hnf_addrow(ZZ_mat<ZT> &B, const vector<Z_NR<ZT>> &w)
{
  /* resize the matrix */
  int r = B.get_rows(), c = B.get_cols();
  B.resize(r + 1, c);

  /* set the initial new row values */
  for (int i = 0; i < c; i++)
  {
    B[r][i] = w[i];
  }

  /* matrix indexes (k = pivot)*/
  r++;
  int i, j, j2, k;
  int l = (c - r) * (c > r);
  /* ZT integers for operations */
  Z_NR<ZT> r1d, r2d, b, u, v, d, q, tmp;

  /* clang-format off */
  _DIAGONAL_INDEX_LOOP_(B,
    /* pre-check : reduce the upper vector with the pivot one */
    if (k) {
      _REDUCE_VECT_XGCD(B,k - 1,j,j2,k,b,d,u,v,r1d,r2d)
    }
    ,
    /* neg : change sign of the row vector if the diagonal entry is negative */
    for (j2 = j; j2 >= 0; j2--) { B[k][j2].neg(B[k][j2]); }
    _REDUCE_LOWER_(B, i, j, j2, k, r);
    ,
    /* zero : modify pivot position */
    { k++; if (l > 0) { l--; } }
    ,
    /* pos : reduce lower entries of column j with row k */
    _REDUCE_LOWER_(B, i, j, j2, k, r);
  )
  /* clang-format on */
}

template <class ZT> int hnf_autoselect(ZZ_mat<ZT> &B)
{
  // rest in development
  // cerr << "warning : xgcd is not suitable for high determinant or big
  // matrices\n";
  return hnf_classical_reduction(B);
}

template <class ZT> int hnf_reduction(ZZ_mat<ZT> &B, HNFMethod method, Z_NR<ZT> &det)
{
  // construct the reduction object
  HNFReduction<ZT> hnf_obj(B, method, det);

  // apply reduction
  hnf_obj.hnf();

  /* Insert any command you desire here after the reduction (for tests or else) */
  // vector<Z_NR<ZT>> v(B.get_cols());
  // for (size_t i = 0; i < v.size(); i++)
  // {
  //   v[i] = rand();
  // }
  // v[0] = 4;
  // v[1] = 1;
  // v[2] = 2;
  // hnf_addrow(B, v);

  // return the status
  return hnf_obj.status;
}

/** enforce instantiation of complete templates **/

// template class HNFReduction<mpz_t>;
template int hnf_reduction(ZZ_mat<mpz_t> &B, HNFMethod method, Z_NR<mpz_t> &det);

FPLLL_END_NAMESPACE
