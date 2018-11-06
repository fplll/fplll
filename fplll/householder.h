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

#ifndef FPLLL_HOUSEHOLDER_H
#define FPLLL_HOUSEHOLDER_H

#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

enum MatHouseholderFlags
{
  HOUSEHOLDER_DEFAULT       = 0,
  HOUSEHOLDER_ROW_EXPO      = 1,
  HOUSEHOLDER_BF            = 2,
  HOUSEHOLDER_OP_FORCE_LONG = 4
};

/**
 * MatHouseholder provides an interface for performing elementary operations on a basis
 * and computing its R matrix using Householder transformations.
 */
template <class ZT, class FT> class MatHouseholder
{
public:
  /**
   * Constructor.
   * The precision of FT must be defined before creating an instance of the
   * class and must remain the same until the object is destroyed (or no longer
   * needed).
   * @param b
   *   The matrix on which row operations are performed. It must not be empty.
   * @param enable_row_expo
   *   If true, each row of b is normalized by a power of 2 before doing
   *   conversion to floating-point, which hopefully avoids some overflows.
   *   This option cannot be enabled if enable_int_gram=true and works only
   *   with FT=double and FT=long double. It is useless and MUST NOT be used
   *   for FT=dpe or FT=mpfr_t.
   * @param enable_bf
   *   Dot product on the basis vector are done on a floatting point version of
   *   it instead of the integer one. bf is refreshed from b when  refresh_R_bf(int)
   *   or refresh_R(int) are called.
   * @param enable_transform
   *   Compute u
   * @param u
   *   If u is not empty, operations on b are also done on u
   *   (in this case both must have the same number of rows).
   *   If u is initially the identity matrix, multiplying transform by the
   *   initial basis gives the current basis.
   * @param enable_inverse_transform
   *   Compute u_inv_t
   * @param u_inv_t
   *   Inverse transform (should be empty, which disables the computation, or
   *   initialized with identity matrix). It works only if u is not empty.
   * @param enable_op_force_long
   *   Affects the behaviour of row_addmul(_we).
   *   See the documentation of row_addmul.
   */
  MatHouseholder(Matrix<ZT> &arg_b, Matrix<ZT> &arg_u, Matrix<ZT> &arg_uinv_t, int flags)
      : b(arg_b), enable_row_expo(flags & HOUSEHOLDER_ROW_EXPO), enable_bf(flags & HOUSEHOLDER_BF),
        enable_transform(arg_u.get_rows() > 0), u(arg_u),
        enable_inverse_transform(arg_uinv_t.get_rows() > 0), u_inv_t(arg_uinv_t),
        row_op_force_long(flags & HOUSEHOLDER_OP_FORCE_LONG)
  {
    // Get the dimensions of the lattice
    d = b.get_rows();
    n = b.get_cols();

    // Any row are known
    n_known_rows = 0;
    n_known_cols = 0;

    // Initialize sigma and V used to compute R
    sigma.resize(d);
    R.resize(d, n);
    // TODO: V does not need to be a matrix, since V is a list of vector of different length
    V.resize(d, n);

    // If enable_bf, initialize bf
    if (enable_bf)
      bf.resize(d, n);

    // Initialize row_expo.
    row_expo.resize(d);
    fill(row_expo.begin(), row_expo.end(), 0);
    // Initialize row_size
    init_row_size.resize(d);
    for (int i = 0; i < d; i++)
      // Capture the shape of b
      init_row_size[i] = max(b[i].size_nz(), 1);

    // Initialize R_history and update_R
    R_history.resize(d);
    for (int i = 0; i < d; i++)
    {
      R_history[i].resize(n);
      for (int j = 0; j < n; j++)
        R_history[i][j].resize(n);
    }
    updated_R = false;

    // Initialize norm_square_b
    norm_square_b.resize(d);
    expo_norm_square_b.resize(d);
    /* fill row_expo with 0, since if enable_row_expo, it will be filled by real value, and
     * otherwise, we essentially 0 - 0 */
    fill(expo_norm_square_b.begin(), expo_norm_square_b.end(), 0);

#ifdef HOUSEHOLDER_PRECOMPUTE_INVERSE
    // Initialize R_inverse_diag
    R_inverse_diag.resize(d);
#endif  // HOUSEHOLDER_PRECOMPUTE_INVERSE

    if (enable_row_expo)
      // Initialize tmp_col_expo
      tmp_col_expo.resize(n);

    /* Initialize values for naively part of the computation
     * Used in is_hlll_reduced at least */
    n_known_rows_naively = 0;
    sigma_naively.resize(d);
    R_naively.resize(d, n);
    V_naively.resize(d, n);
    row_expo_naively.resize(d);
    fill(row_expo_naively.begin(), row_expo_naively.end(), 0);

#ifdef DEBUG
    for (int i = 0; i < d; i++)
    {
      for (int j = 0; j < n; j++)
      {
        V(i, j).set_nan();
        V_naively(i, j).set_nan();
      }
    }
#endif  // DEBUG
  }

  ~MatHouseholder() {}

  /**
   * Returns f (* 2^expo if enable_row_expo) = R(i, j).
   *
   * Returns reference to `f`.
   */
  inline void get_R(FT &f, int i, int j, long &expo);

  /**
   * Sets R(i, j) to f.
   */
  inline void set_R(FT &f, int i, int j);

  /**
   * Returns R[i].
   */
  inline MatrixRow<FT> get_R(int i, long &expo);

  /**
   * Returns the R matrix
   * expo is set to row_expo
   */
  const Matrix<FT> &get_R(vector<long> &expo)
  {
    expo = row_expo;
    return R;
  }

  /**
   * Returns b[i].
   */
  MatrixRow<ZT> get_b(int i);

  /**
   * Returns the b matrix
   */
  const Matrix<ZT> &get_b() { return b; }

  /**
   * Apply Householder transformation on row i for columns [0, i).
   * If last_j, apply Householder transformation on row i, from cols [0, i].
   */
  void update_R(int i, bool last_j);

  /**
   * Apply Householder transformation on row i.
   */
  void update_R(int i);

  /**
   * Finalize Householder transformation on row i (especially after update_R(i, false))
   */
  void update_R_last(int i);

  /**
   * Full computation of the matrix R.
   */
  inline void update_R();

  /**
   * Retun the dimensions of the lattice
   */
  inline int get_d() { return d; }
  inline int get_n() { return n; }

  /**
   * Compute the squared norm of b[k].
   */
  inline void norm_square_b_row(FT &f, int k, long &expo);

  /**
   * Truncated squared norm of R[k], with coefficients of R[k][beg..end-1].
   */
  inline void norm_square_R_row(FT &f, int k, int beg, int end, long &expo);

  // b[k] = b[k] + ZT(xf) * b[i]
  // R[k] = R[k] + xf * R[i]
  // Not necessary to make transformation on bf, since basically, after each size_reduce, we call
  // refresh_R_bf, which set bf to the correct value from b directly.
  // xf must be non-zero.
  void size_reduce(const FT &xf, int k, int i);

  /**
   * Swap row i and j of b, bf, R, V, u and u_inv_t
   * Swap element i and j in sigma, row_expo, norm_square_b, expo_norm_square_b, init_row_size and
   * R_history.
   */
  void swap(int i, int j);

  /**
   * Update n_known_rows to k.
   */
  inline void invalidate_row(int k);

  /**
   * Return values enable_row_expo and enable_bf
   */
  inline bool is_enable_row_expo() { return enable_row_expo; }
  inline bool is_enable_bf() { return enable_bf; }

  /**
   * Return value of updated_R
   */
  inline bool get_updated_R() { return updated_R; }
  /**
   * Set updated_R to false
   * updated_R is set to true in recover_R
   */
  inline void set_updated_R_false() { updated_R = false; }

  /**
   * Get the precomputation of R(i, i)
   */
  inline FT get_R_inverse_diag(int i) { return R_inverse_diag[i]; }

  /*
   * Recover R[i] from the precomputed values of R stored in R_history, i.e., R[i] is correct for
   * the coefficients [0,
   * i) and a call to update_R_last(i) will fully compute R[i].
   */
  inline void recover_R(int i);

  /**
   * Return row_expo[i]
   */
  inline long get_row_expo(int i) { return row_expo[i]; }

  /**
   * Returns the value of row_op_force_long
   */
  inline bool is_row_op_force_long() { return row_op_force_long; }

  /**
   * Set bf[i] and R[i] to b[i].
   * Precompute square norm of b[i].
   */
  void refresh_R_bf(int i);
  /**
   * Set bf and R to b.
   * Precompute squared norm of all the vectors of b
   */
  inline void refresh_R_bf();

  /**
   * Set R[i] to b[i].
   */
  void refresh_R(int i);
  /**
   * Set R to b.
   */
  inline void refresh_R();

  /**
   * Set in f the precomputed squared norm of b[i] and in expo the exponent such that ||b[i]||^2 = f
   * * 2^{expo}
   */
  inline void get_norm_square_b(FT &f, int i, long &expo);

private:
  /**
   * Number of rows of b (dimension of the lattice).
   */
  int d;

  /**
   * Number of columns of b (dimension of the lattice).
   */
  int n;

  /**
   * Basis of the lattice
   */
  Matrix<ZT> &b;

  /**
   * b = R * q_householder (q_householder is not computed)
   * R is lower triangular and the diagonal coefficient are >= 0.
   */
  Matrix<FT> R;

  /**
   * Vector v following [MSV, ISSAC'09].
   */
  Matrix<FT> V;

  /**
   * Sigma values following [MSV, ISSAC'09].
   */
  vector<FT> sigma;

  /**
   * R[i] is invalid for i >= n_known_rows.
   */
  int n_known_rows;

  /** Normalization of each row of b by a power of 2. */
  const bool enable_row_expo;

  /**
   * When enable_row_expo=true, row_expo[i] is the smallest non-negative integer
   * such that b(i, j) &lt;= 2^row_expo[i] for all j. Otherwise this array is empty.
   */
  vector<long> row_expo;

  /* Used by update_R. Temporary variable. */
  vector<long> tmp_col_expo;

  // Temporary variables
  FT ftmp0, ftmp1, ftmp2, ftmp3;
  ZT ztmp0, ztmp1;

  // init_row_size[i] = (last non-zero column in the i-th row of b) + 1
  vector<int> init_row_size;
  // n_known_cols (last non-zero column of the discovered rows) + 1
  int n_known_cols;

  /**
   * b[i] := b[i] + x * 2^expo_add * b[j].
   * Special cases |x| &lt;= 1 and |x| &lt;= LONG_MAX are optimized.
   * x should be a non-zero FT. x*2^expo_add must represent an integer (typically coming from
   * rnd_we).
   * If row_op_force_long=true, x is always converted to (2^expo * long) instead
   * of (2^expo * ZT), which is faster if ZT=mpz_t but might lead to a loss of
   * precision.
   */
  void row_addmul_we(int i, int j, const FT &x, long expo_add);
  /**
   * Special cases of row_addmul_we
   */
  void row_add(int i, int j);
  void row_sub(int i, int j);
  void row_addmul_si(int i, int j, long x);
  void row_addmul_si_2exp(int i, int j, long x, long expo);
  void row_addmul_2exp(int i, int j, const ZT &x, long expo);

  /**
   * Basis of the lattice (floatting point)
   */
  Matrix<FT> bf;

  /**
   * Do we compute dot_product of bf instead of b?
   */
  const bool enable_bf;

  /**
   * R_history stores the history of the computation of R
   * Example: R[i][j][k] is the index k of R[i] when the coefficient j was known
   */
  vector<vector<vector<FT>>> R_history;

  // If updated_R, R[n_known_rows][0] to R[n_known_rows][n_known_rows - 1] is valid
  bool updated_R;

  // Stores at index i the inverse of R(i, i)
  vector<FT> R_inverse_diag;

  // Compute the unimodular matrix u
  const bool enable_transform;

  Matrix<ZT> &u;  // Transform

  /**
   * Computation of the inverse transform matrix (transposed).
   * This works only if enable_transform=true.
   * This matrix has very large coefficients, computing it is slow.
   */
  const bool enable_inverse_transform;

  Matrix<ZT> &u_inv_t;  // Transposed inverse transform

  /**
   * Changes the behaviour of row_addmul(_we).
   * See the description of row_addmul.
   */
  const bool row_op_force_long;

  // Store the approximate norm of b[i].
  vector<FT> norm_square_b;
  // If enable_row_expo, ||b[i]||^2 = norm_square_b[i] * 2^{expo_norm_square_b}
  vector<long> expo_norm_square_b;

  /* Objects and methods for the naive computation of the R factor using Householder. */

public:
  /**
   * Full computation of the matrix R.
   */
  void update_R_naively();

  /**
   * Apply Householder transformation on row i.
   */
  void update_R_naively(int i);

  /**
   * Return R_naively(i, j) = f (* 2^expo, if enable_row_expo)
   */
  inline void get_R_naively(FT &f, int i, int j, long &expo);

  /**
   * Squared norm of b[k].
   * Use row_expo_naively if enable_row_expo is used.
   * f * 2^expo = ||b[i]||^2
   */
  inline void norm_square_b_row_naively(FT &f, int k, long &expo);

  /**
   * Truncated norm square of R_naively[k], with coefficients of R_naively[k][0..end-1].
   */
  inline void norm_square_R_row_naively(FT &f, int k, int end, long &expo);

private:
  /**
   * b = R * q_householder.
   * R is lower triangular and the diagonal coefficient are >= 0.
   */
  Matrix<FT> R_naively;

  /**
   * Vector v following [MSV, ISSAC'09].
   */
  Matrix<FT> V_naively;

  /**
   * Sigma values following [MSV, ISSAC'09].
   */
  vector<FT> sigma_naively;

  /**
   * When enable_row_expo=true, row_expo_naively[i] is the smallest non-negative integer
   * such that b(i, j) &lt;= 2^row_expo_naively[i] for all j. Otherwise this array is empty.
   */
  vector<long> row_expo_naively;

  /**
   * R[i] is invalid for i >= n_known_rows_naively.
   */
  int n_known_rows_naively;
};

template <class ZT, class FT>
inline void MatHouseholder<ZT, FT>::get_R(FT &f, int i, int j, long &expo)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d && j >= 0 && j <= i);
  f    = R(i, j);
  expo = row_expo[i];
}

template <class ZT, class FT> inline MatrixRow<FT> MatHouseholder<ZT, FT>::get_R(int i, long &expo)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d);
  expo = row_expo[i];

  return R[i];
}

template <class ZT, class FT> MatrixRow<ZT> MatHouseholder<ZT, FT>::get_b(int i)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d);
  return b[i];
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::set_R(FT &f, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d && i >= j && j >= 0);
  FPLLL_DEBUG_CHECK(j <= i);
  R(i, j) = f;
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::update_R(int i)
{
  // last_j = true since we want the full computation of R[i]
  update_R(i, true);
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::update_R()
{
  for (int i = 0; i < d; i++)
    update_R(i);
}

template <class ZT, class FT>
inline void MatHouseholder<ZT, FT>::norm_square_b_row(FT &f, int k, long &expo)
{
  FPLLL_DEBUG_CHECK(k >= 0 && k < d);

  if (enable_bf)
  {
    bf[k].dot_product(f, bf[k], 0, n_known_cols);

    if (enable_row_expo)
      expo = 2 * row_expo[k];
    else
      expo = 0;
  }
  else
  {
    b[k].dot_product(ztmp0, b[k], 0, n_known_cols);

    if (enable_row_expo)
      ztmp0.get_f_exp(f, expo);
    else
    {
      f.set_z(ztmp0);
      expo = 0;
    }
  }
}

template <class ZT, class FT>
inline void MatHouseholder<ZT, FT>::norm_square_R_row(FT &f, int k, int beg, int end, long &expo)
{
  FPLLL_DEBUG_CHECK(k >= 0 && k < d);
  FPLLL_DEBUG_CHECK(beg <= end && end <= k);
  if (end == beg)
    f = 0.0;
  else
    R[k].dot_product(f, R[k], beg, end);

  if (enable_row_expo)
    expo = 2 * row_expo[k];
  else
    expo = 0;
}

// TODO: test seems to be strange
template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::invalidate_row(int k)
{
  if (k < n_known_rows)
    n_known_rows = k;
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::recover_R(int i)
{
  for (int k = 0; k < i - 1; k++)
    R(i, k) = R_history[i][k][k];
  for (int k = i - 1; k < n; k++)
    R(i, k) = R_history[i][i - 1][k];

  updated_R = true;
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::refresh_R_bf()
{
  for (int i = 0; i < d; i++)
    refresh_R_bf(i);
}

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::refresh_R()
{
  for (int i = 0; i < d; i++)
    refresh_R(i);
}

template <class ZT, class FT>
inline void MatHouseholder<ZT, FT>::get_norm_square_b(FT &f, int i, long &expo)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d);
  expo = expo_norm_square_b[i];
  f    = norm_square_b[i];
}

/* Objects and methods for the naive computation of the R factor using Householder. */

template <class ZT, class FT> inline void MatHouseholder<ZT, FT>::update_R_naively()
{
  for (int i = 0; i < d; i++)
    update_R_naively(i);
}

/* TODO: refactorize. */
template <class ZT, class FT>
inline void MatHouseholder<ZT, FT>::get_R_naively(FT &f, int i, int j, long &expo)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < d && j >= 0 && j <= i);
  f    = R_naively(i, j);
  expo = row_expo_naively[i];
}

template <class ZT, class FT>
inline void MatHouseholder<ZT, FT>::norm_square_R_row_naively(FT &f, int k, int end, long &expo)
{
  FPLLL_DEBUG_CHECK(k >= 0 && k < d);
  FPLLL_DEBUG_CHECK(0 <= end && end <= k);
  if (end == 0)
    f = 0.0;
  else
    R_naively[k].dot_product(f, R_naively[k], 0, end);

  if (enable_row_expo)
    expo = 2 * row_expo_naively[k];
  else
    expo = 0;
}

template <class ZT, class FT>
inline void MatHouseholder<ZT, FT>::norm_square_b_row_naively(FT &f, int k, long &expo)
{
  FPLLL_DEBUG_CHECK(k >= 0 && k < d);
  if (enable_row_expo)
  {
    b[k].dot_product(ztmp0, b[k], 0, n);
    ztmp0.get_f_exp(f, expo);
  }
  else
  {
    expo = 0;
    b[k].dot_product(ztmp0, b[k], 0, n);
    f.set_z(ztmp0);
  }
}

FPLLL_END_NAMESPACE

#endif  // FPLLL_HOUSEHOLDER_H
