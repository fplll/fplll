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

// Class method templates

bool HNFReduction::set_status(int new_status)
{
  status = new_status;
  return status == RED_SUCCESS;
}

void HNFReduction::is_reduced()
{
  set_status(is_hnf_reduced(b));
}

void HNFReduction::hnf()
{
  int new_status = RED_HNF_FAILURE;
  if (method == HM_AUTO)
  {
    new_status = hnf_autoselect(b);
  }
  else if (method == HM_XGCD)
  {
    new_status = hnf_xgcd_reduction(b);
  }
  else if (method == HM_CLASSIC)
  {
    new_status = hnf_classical_reduction(b);
  }
  else if (method == HM_MINORS)
  {
    new_status = hnf_minors_reduction(b);
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
  //is_reduced();
  set_status(new_status);
}

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
        if (B[i][j] != 0) { return RED_HNF_FAILURE;}                                                             \
      }                                                                                            \
      ,                                                                                            \
      /* neg : diagonal shouldn't be negative */                                                   \
      { return RED_HNF_FAILURE; }                                                                                \
      ,                                                                                            \
      /* zero : pivot row position doesn't change, increment to compensate */                      \
      /* lower the limit as well, we skipped one column without increasing row */                  \
      { k++; if (l > 0) { l--; } }                                                                 \
      ,                                                                                            \
      /* pos : checks if everything below is also positive and lower than the pivot */             \
      /* leave a potential extra instruction */                                                    \
      for (i = k + 1; i < r; i++) {                                                                \
        if ((B[i][j] > B[k][j]) || (B[i][j] < 0)) { return RED_HNF_FAILURE; }                                    \
      }                                                                                            \
      /* leave a potential extra instruction */                                                    \
      {INSTRUCTION_AT_PIVOT}                                                                       \
    )                                                                                              \
  }

/* clang-format on */

// HNF computation macro, reduce lower coefficients of column j by the pivot in row k
#define _REDUCE_LOWER_(B, i, j, j2, k, r)                                                          \
  {                                                                                                \
    for (i = k + 1; i < r; i++)                                                                    \
    {                                                                                              \
      q.fdiv_q(B[i][j], B[k][j]);                                                                  \
      for (j2 = j; j2 >= 0; j2--)                                                                  \
      {                                                                                            \
        B[i][j2].submul(q, B[k][j2]);                                                              \
      }                                                                                            \
    }                                                                                              \
  }

// Takes a matrix and two row indexes and applies the xgcd reduction, and updates pivot
// j is the head coefficient's column. d,u,v,r1d,r2d are variables used for xgcd
#define _REDUCE_VECT_XGCD(B, to_red, j, j2, pivot,b,d,u,v,r1d,r2d)                                                          \
  {                                                                                                \
    /* skip zero */                                                                                \
    if (B[to_red][j] == 0) { continue; }                                                                \
    /* reduce row to_red with row pivot */                                                                  \
    d.xgcd(u, v, B[pivot][j], B[to_red][j]);                                                                \
    /* This is in case the pivot is the gcd, saves some time (useful for ones for sure) */\
    if (mpz_cmpabs(d.get_data(), B[pivot][j].get_data()) == 0)\
    {\
      b.divexact(B[to_red][j], B[pivot][j]);\
      for (j2 = j; j2 >= 0; j2--)\
      {\
        B[to_red][j2].submul(b, B[pivot][j2]);\
      }\
      continue;\
    }\
    r2d.divexact(B[to_red][j], d);                                                                      \
    r1d.divexact(B[pivot][j], d);                                                                      \
    /* only compute relevant values (rest should be guaranteed zero) */                            \
    for (j2 = j; j2 >= 0; j2--) {                                                                  \
      /* precompute values */                                                                      \
      b.mul(u, B[pivot][j2]);                                                                          \
      b.addmul(v, B[to_red][j2]);                                                                       \
      /* new vector i-1 value */                                                                   \
      B[to_red][j2].mul(r1d, B[to_red][j2]);                                                                 \
      B[to_red][j2].submul(r2d, B[pivot][j2]);                                                              \
      /* new vector i value */                                                                     \
      B[pivot][j2] = b;                                                                                \
    }                                                                                              \
  }

int is_hnf_reduced(const ZZ_mat<mpz_t> &B)
{
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols(); /* matrix indexes (k = pivot)*/
  int i, j, k;                            /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);              /* main loop, decreases from last column to the limit */
  _HNF_CHECK_(B, {});
  return RED_SUCCESS;
}

int in_lattice_given_hnf(ZZ_mat<mpz_t> &B, const vector<Z_NR<mpz_t>> &w)
{
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols();
  /* matrix indexes (k = pivot)*/
  int i, j, j2, k;
  /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);

  vector<Z_NR<mpz_t>> v = w;
  Z_NR<mpz_t> q;

  /* test if the vector size is correct */
  if ((int)v.size() != r)
  {
    cerr << "in_hnf error : matrix-vector sizes do not match\n";
    return -1;
  }
  /* clang-format off */
  _HNF_CHECK_(B,
    /* reduces the vector by the row pivot vector for membership test */
    q.fdiv_q(v[j], B[k][j]);
    for (j2 = j; j2 >= 0; j2--)
      {
        v[j2].submul(q, B[k][j2]);
      }
  );
  /* clang-format on */

  // membership test : the vector after reduction should be a zero one.
  for (j = 0; j < c; j++)
  {
    if (v[j] != 0)
    {
      return -1;
    }
  }

  return 0;
}

int in_lattice_given_hnf(ZZ_mat<mpz_t> &B, const ZZ_mat<mpz_t> &A)
{

  ZZ_mat<mpz_t> tmp = A;
  Z_NR<mpz_t> q;

  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols(); /* matrix indexes (k = pivot)*/
  int i, j, j2, k;                        /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);              /* main loop, decreases from last column to the limit */

  /* test if the matrix size is correct */
  if ((A.get_cols() != B.get_cols()) || (A.get_rows() != B.get_rows()))
  {
    cerr << "in_hnf error : matrix-matrix sizes do not match\n";
    return -1;
  }

  /* clang-format off */
  _HNF_CHECK_(B,
    /* reduces the matrix tmp by the row pivot vector of A for membership test */
    for (i = 0; i < r; i++)
      {
       q.fdiv_q(tmp[i][j], B[k][j]);
       for (j2 = j; j2 >= 0; j2--)
         {
           tmp[i][j2].submul(q, B[k][j2]);
         }
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
        return -1;
      }
    }
  }
  return 0;
}

int hnf_xgcd_reduction(ZZ_mat<mpz_t> &B)
{
  /* matrix indexes (k = pivot)*/
  int i, j, j2, k;
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols();
  /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);
  /* ZT integers for operations */
  Z_NR<mpz_t> r1d, r2d, b, u, v, d, q;

  /* initializes the transformation matrix */
  // ZZ_mat<mpz_t> U;
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

int hnf_classical_reduction(ZZ_mat<mpz_t> &B)
{
  /* matrix indexes (k = pivot)*/
  int i, j, j2, k, i_min;
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols();
  /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r);

  /* ZT integers for operations */
  Z_NR<mpz_t> q;
  /* value to check whether or not a column is already reduced */
  int reduced;

  /* initializes the transformation matrix */
  // ZZ_mat<mpz_t> U;
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
      if ( B[k][j] != 0 && mpz_cmpabs(B[i][j].get_data(), B[k][j].get_data()) > 0)
        {
          i_min = k;
        } else {
          i_min = i;
        }
      /* have to find the row with the minimum non-zero absolute value */
      for (i = i - 1 ; i >= 0; i--) {
        if (B[i][j] != 0 && mpz_cmpabs(B[i_min][j].get_data(), B[i][j].get_data()) > 0)
          { i_min = i; }
      }
      /* use the minimum row as pivot, reduce the rest above */
      B.swap_rows(i_min, k);
      for (i = k - 1, reduced = 1; (i >= 0); i--)
      {
        q.fdiv_q(B[i][j], B[k][j]);
        for (j2 = j; j2 >= 0; j2--)
        {
          B[i][j2].submul(q, B[k][j2]);
        }
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

int hnf_modular_reduction(ZZ_mat<mpz_t> &B, const Z_NR<mpz_t> D)
{
  /* matrix indexes (k = pivot)*/
  int i, j, k;
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols();
  if (r < c)
  {
    cerr << "modular method requires at least as many rows as columns" << '\n';
    return -1;
  }

  /* ZT integers for operations */
  Z_NR<mpz_t> R = D, R2, d, u, v, r1d, r2d, b, q;

  /* initializes the transformation matrix */
  // ZZ_mat<mpz_t> U;
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
      q.fdiv_q(B[i][k], B[k][k]);
      for (j = k; j >= 0; j--)
      {
        B[i][j].submul(q, B[k][j]);
      }
    }

    /* Sets the new determinant remainder */
    R.divexact(R, d);
  }

  // return in_lattice_given_hnf(B,U);
  return 0;
}

int hnf_minors_reduction(ZZ_mat<mpz_t> &B)
{
  /* matrix indexes (p_i, p_j = pivot being reduced, k = block_size)*/
  int k, j, j2, l, block_size, max_minor, p_i, p_j;
  /* matrix bounds */
  int r = B.get_rows(), c = B.get_cols();
  /* ZT integers for operations */
  if (r < c)
  {
    max_minor = r;
  } else {
    max_minor = c;
  }
  /* this variable indicates the row to swap with the pivot */
  int i_swap = 0;

  Z_NR<mpz_t> r1d, r2d, b, u, v, d, q;

  /* initializes the transformation matrix */
  // ZZ_mat<mpz_t> U;
  // U = B;
  // U.gen_identity(B.get_rows());

  /* put the kth principal minor in HNF (starting from bottom right corner)
   * but without computing the lower triangle for now (get diagonal coefficients)
   * k indicates the block size */
  for (block_size = 1; block_size <= max_minor ; block_size++)
  {
    p_i = r - block_size;
    p_j = c - block_size;
    l = ((c - block_size) * (c > block_size)) + 1;
    // for (j = (c - 1), k = (r - 1); j > (c - block_size) && k > (r - block_size) ; j--, k--)
    for (j = (c - 1), k = (r - 1); j >=l ; j--, k--)
    {
      /* reduces p_i with row i */
      _REDUCE_VECT_XGCD(B, p_i, j, j2, k,b,d,u,v,r1d,r2d)
    }
    /* if the next last pivot is zero we swap row k for some other row (starting with the first) */
    if (B[p_i][p_j] == 0)
    {
      if (i_swap != p_i) {
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
  // for (j = (c - 1), k = (r - 1); j >= c - max_minor ; j--, k--)
  for (j = (c - 1), k = (r - 1); j >=l ; j--, k--)
  {
    /* reduce extra rows */
    for (block_size = 0; block_size < r - max_minor; block_size++)
    {
      _REDUCE_VECT_XGCD(B, block_size, j, j2, k,b,d,u,v,r1d,r2d)
    }

    /* reduce below diagonal elements */
    _REDUCE_LOWER_(B,p_i,j,j2,k,r);
  }

  return 0;
}

int hnf_reduction(ZZ_mat<mpz_t> &B, HNFMethod method, Z_NR<mpz_t> &det)
{
  //construct the reduction object
  HNFReduction hnf_obj(B,method,det);

  //apply reduction
  hnf_obj.hnf();

  //return the status
  return hnf_obj.status;
}

int hnf_autoselect(ZZ_mat<mpz_t> &B)
{
  // rest in development
  // cerr << "warning : xgcd is not suitable for high determinant or big
  // matrices\n";
  return hnf_classical_reduction(B);
}

FPLLL_END_NAMESPACE
