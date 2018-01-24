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

// HNF checking macro to check if B is a HNF, and add an instruction if the pivot is correct
#define _HNF_CHECK_(B, INSTRUCTION_AT_PIVOT)                                                       \
  {                                                                                                \
    _DIAGONAL_INDEX_LOOP_(B, /* checks if everything above diagonal is zero */                     \
                          for (i = k - 1; i >= 0; i--) {                                           \
                            if (B[i][j] != 0)                                                      \
                            {                                                                      \
                              return 1;                                                            \
                            }                                                                      \
                          }, /* neg : diagonal shouldn't be negative */                            \
                          {                                                                        \
                            return 1;                                                              \
                          }, /* zero : pivot row position doesn't change, increment to compensate  \
                                */ /* lower the limit as well, we skipped one column without       \
                                      increasing row */                                            \
                          {                                                                        \
                            k++;                                                                   \
                            if (l > 0)                                                             \
                            {                                                                      \
                              l--;                                                                 \
                            }                                                                      \
                          }, /* pos : checks if everything below is also positive and lower than   \
                                the pivot */ /* leave a potential extra instruction */             \
                          for (i = k + 1; i < r; i++) {                                            \
                            if ((B[i][j] > B[k][j]) || (B[i][j] < 0))                              \
                            {                                                                      \
                              return 1;                                                            \
                            }                                                                      \
                          } {INSTRUCTION_AT_PIVOT})                                                \
  }

// HNF computation macro, reduce lower coefficients of column j by the pivot in row k
#define _REDUCE_LOWER_(B, i, j2)                                                                   \
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

FPLLL_BEGIN_NAMESPACE

int is_hnf_reduced(ZZ_mat<mpz_t> &B)
{
  /* matrix bounds */
  int r = B.get_rows(),
      c = B.get_cols();      /* matrix indexes (k = pivot)*/
  int i, j, k;               /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r); /* main loop, decreases from last column to the limit */
  _HNF_CHECK_(B, {});
  return 0;
}

int in_lattice_given_hnf(ZZ_mat<mpz_t> &B, const vector<Z_NR<mpz_t>> &w)
{
  /* matrix bounds */
  int r = B.get_rows(),
      c = B.get_cols();      /* matrix indexes (k = pivot)*/
  int i, j, j2, k;           /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r); /* main loop, decreases from last column to the limit */

  vector<Z_NR<mpz_t>> v = w;
  Z_NR<mpz_t> q;

  /* test if the vector size is correct */
  if ((int)v.size() != r)
  {
    cerr << "in_hnf error : matrix-vector sizes do not match\n";
    return -1;
  }

  _HNF_CHECK_(B, { /* reduces the vector by the row pivot vector for membership test */
                   q.fdiv_q(v[j], B[k][j]);
                   for (j2 = j; j2 >= 0; j2--)
                   {
                     v[j2].submul(q, B[k][j2]);
                   }
  });

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
  int r = B.get_rows(),
      c = B.get_cols();      /* matrix indexes (k = pivot)*/
  int i, j, j2, k;           /* the limit "l" depends of matrix dimensions */
  int l = (c - r) * (c > r); /* main loop, decreases from last column to the limit */

  /* test if the matrix size is correct */
  if ((A.get_cols() != B.get_cols()) || (A.get_rows() != B.get_rows()))
  {
    cerr << "in_hnf error : matrix-matrix sizes do not match\n";
    return -1;
  }

  _HNF_CHECK_(B, { /* reduces the matrix tmp by the row pivot vector of A for membership test */
                   for (i = 0; i < r; i++)
                   {
                     q.fdiv_q(tmp[i][j], B[k][j]);
                     for (j2 = j; j2 >= 0; j2--)
                     {
                       tmp[i][j2].submul(q, B[k][j2]);
                     }
                   }
                   cout << "HNF reduction for pivot " << k << endl << B << endl;
                   cout << "tmp reduction for pivot " << k << endl << tmp << endl;
  });

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

  _DIAGONAL_INDEX_LOOP_(
      B, /* pre-check : iterates above the pivot to construct it */
      for (i = k - 1; i >= 0;
           i--) { /* skip zeroes (if any) */
                  if (B[i + 1][j] == 0)
                  {
                    continue;
                  } /* reduce row i + 1 with row i */
                  d.xgcd(u, v, B[i][j], B[i + 1][j]);
                  r2d.divexact(B[i + 1][j], d);
                  r1d.divexact(
                      B[i][j],
                      d); /* only compute relevant values (rest should be guaranteed zero) */
                  for (j2 = j; j2 >= 0; j2--)
                  { /* precompute values */
                    b.mul(u, B[i][j2]);
                    b.addmul(v, B[i + 1][j2]); /* new vector i-1 value */
                    B[i + 1][j2].mul(r1d, B[i + 1][j2]);
                    B[i + 1][j2].submul(r2d, B[i][j2]); /* new vector i value */
                    B[i][j2] = b;
                  }
      } /* swap first row with the pivot row */
      B.swap_rows(0, k);
      , /* neg : change sign of the row vector if the diagonal entry is negative */
      for (j2 = j; j2 >= 0; j2--) { B[k][j2].neg(B[k][j2]); }, /* zero : modify pivot position */
      {
        k++;
        if (l > 0)
        {
          l--;
        }
      }, /* pos : reduce lower entries of column j with row k */
      _REDUCE_LOWER_(B, i, j2))

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
  Z_NR<mpz_t> min, q;
  /* value to check whether or not a column is already reduced */
  int reduced;

  /* initializes the transformation matrix */
  // ZZ_mat<mpz_t> U;
  // U = B;
  // U.gen_identity(B.get_rows());

  _DIAGONAL_INDEX_LOOP_(
      B, /* pre-check : iterates above the pivot to construct it */ /* check if the column j is
                                                                       already reduced above*/
      for (i = k - 1, reduced = 1; (i >= 0) && reduced; i--) {
        reduced = (B[i][j] == 0);
      } if (!reduced) { /* have to find the row with the minimum non-zero absolute value*/
                        for (i_min = k, i = k - 1; i >= 0; i--)
                        {
                          if (mpz_cmpabs(B[i_min][j].get_data(), B[i][j].get_data()) > 0 &&
                              B[i][j] != 0)
                          {
                            i_min = i;
                          }
                        } /* use the minimum row as pivot, make it positive and reduce the rest */
                        B.swap_rows(i_min, k);
                        for (i = k - 1, reduced = 1; (i >= 0); i--)
                        {
                          q.fdiv_q(B[i][j], B[k][j]);
                          for (j2 = j; j2 >= 0; j2--)
                          {
                            B[i][j2].submul(q, B[k][j2]);
                          }
                        } /* have to check again if we're finished before going on with the checks*/
                        j++;
                        k++;
                        continue;
      }, /* neg : change sign of the row vector if the diagonal entry is negative, then reduce */
      for (j2 = j; j2 >= 0; j2--) {
        B[k][j2].neg(B[k][j2]);
      } _REDUCE_LOWER_(B, i, j2), /* zero : modify pivot position */
      {
        k++;
        if (l > 0)
        {
          l--;
        }
      }, /* pos : reduce lower entries of column j with row k */
      _REDUCE_LOWER_(B, i, j2))

  // return in_lattice_given_hnf(B,U);
  return 0;
}

int hnf(ZZ_mat<mpz_t> &B, HNFMethod method)
{
  // rest in development
  if (method == HM_XGCD)
  {
    return hnf_xgcd_reduction(B);
  }
  else if (method == HM_CLASSIC)
  {
    return hnf_classical_reduction(B);
  }
  else
  {
    cerr << "HNF method not implemented yet\n";
    return -1;
  }
}

int hnf_autoselect(ZZ_mat<mpz_t> &B)
{
  // rest in development
  // cerr << "warning : xgcd is not suitable for high determinant or big
  // matrices\n";
  return hnf(B, HM_CLASSIC);
}

FPLLL_END_NAMESPACE
