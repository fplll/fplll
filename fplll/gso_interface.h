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

#ifndef FPLLL_GSOInterface_H
#define FPLLL_GSOInterface_H

#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

enum MatGSOInterfaceFlags
{
  GSO_DEFAULT       = 0,
  GSO_INT_GRAM      = 1,
  GSO_ROW_EXPO      = 2,
  GSO_OP_FORCE_LONG = 4,
  GSO_HOUSEHOLDER   = 16  // GSO_GIVENS in kodebro:master_with_givens is equal to 8
};

/**
   @brief Use Gaussian Heuristic to compute a bound on the length of the
   shortest vector

   @param max_dist         output
   @param max_dist_expo    exponent of output
   @param block_size       block size
   @param root_det         root determinant of lattice
   @param gh_factor        factor by which to multiple bound

   @return new bound if `gh_factor * GH` is shorter than `max_dist`, otherwise `max_dist` is
   unchanged.
*/

template <class FT>
void adjust_radius_to_gh_bound(FT &max_dist, long max_dist_expo, int block_size, const FT &root_det,
                               double gh_factor);

/**
 * MatGSOInterface provides an interface for performing elementary operations on a basis
 * and computing its Gram matrix and its Gram-Schmidt orthogonalization.
 * The Gram-Schmidt coefficients are computed on demand. The object keeps track
 * of which coefficients are valid after each row operation.
 */

template <class ZT, class FT> class MatGSOInterface
{
public:
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
  virtual ~MatGSOInterface();

  MatGSOInterface(Matrix<ZT> &arg_u, Matrix<ZT> &arg_uinv_t, int flags)
      : enable_int_gram(flags & GSO_INT_GRAM), enable_row_expo(flags & GSO_ROW_EXPO),
        enable_transform(arg_u.get_rows() > 0), enable_inverse_transform(arg_uinv_t.get_rows() > 0),
        enable_householder(flags & GSO_HOUSEHOLDER), row_op_force_long(flags & GSO_OP_FORCE_LONG),
        u(arg_u), u_inv_t(arg_uinv_t), n_known_rows(0), n_source_rows(0), n_known_cols(0),
        cols_locked(false), alloc_dim(0), gptr(nullptr)
  {
#ifdef DEBUG
    row_op_first = row_op_last = -1;
#endif
  }

  /**
   * Number of rows of b (dimension of the lattice).
   * Can be changed with create_row or remove_last_row.
   */
  int d;

  /**
   * Basis of the lattice
   */
  // Matrix<ZT> &b;

  /**
   * The next five functions make calls
   * from lll.cpp and bkz.cpp indirect.
   */

  /**
   * Returns maximum exponent of b. In the gram
   * version it returns a half times the
   * maximum exponent of g.
  */
  virtual long get_max_exp_of_b() = 0;

  /** Returns true if the ith row
   *  of b is zero. In the gram version
   * it returns true if g(i,i) is zero.
   */
  virtual bool b_row_is_zero(int i) = 0;

  /** Returns number of columns of b. In
   *  the gram version it returns the number
   *  of columns of g.
   */
  virtual int get_cols_of_b() = 0;

  /** Returns number of rows of b. In
   * the gram version it returns the number
   * of of rows of g. This function is made
   * to reduce code repetition (dump_mu/dump_r)
   */
  virtual int get_rows_of_b() = 0;

  /** Negates the ith row of b. Needed
    * by dbkz_postprocessing.
    */
  virtual void negate_row_of_b(int i) = 0;

  /**
   * When enable_row_expo=true, row_expo[i] is the smallest non-negative integer
   * such that b(i, j) &lt;= 2^row_expo[i] for all j. Otherwise this array is empty.
   */
  vector<long> row_expo;

  /**
   * Must be called before a sequence of row_addmul(_we).
   */
  inline void row_op_begin(int first, int last);

  /**
   * Must be called after a sequence of row_addmul(_we). This invalidates the
   * i-th line of the GSO.
   */
  void row_op_end(int first, int last);

  /**
   * Returns Gram matrix coefficients (0 &lt;= i &lt; n_known_rows and
   * 0 &lt;= j &lt;= i).
   * If enable_row_expo=false, returns the dot product (b[i], b[j]).
   * If enable_row_expo=true, returns
   * (b[i], b[j]) / 2 ^ (row_expo[i] + row_expo[j]).
   *
   * Returns reference to `f`.
   */
  virtual FT &get_gram(FT &f, int i, int j) = 0;

  /**
   * Returns the mu matrix
   * Coefficients of the Gram Schmidt Orthogonalization
   * (lower triangular matrix)
   * mu(i, j) = r(i, j) / ||b*_j||^2.
   */
  const Matrix<FT> &get_mu_matrix() { return mu; }

  /**
   * Returns the r matrix
   * Coefficients of the Gram Schmidt Orthogonalization
   * (lower triangular matrix)
   */
  const Matrix<FT> &get_r_matrix() { return r; }

  /**
   * Returns the g matrix (Z_NR version of r)
   * Coefficients of the Gram Schmidt Orthogonalization
   * (lower triangular matrix)
   */
  const Matrix<ZT> &get_g_matrix()
  {
    if (gptr == nullptr)
    {
      throw std::runtime_error("Error: gptr == nullpointer.");
    }
    return *gptr;
  }

  /**
   * Returns f = mu(i, j) and expo such that
   * f * 2^expo = (b_i, b*_j) / ||b*_j||^2.
   * If enable_row_expo=false, expo is always 0.
   * If enable_row_expo=true, expo = row_expo[i] - row_expo[j]
   * It is assumed that mu(i, j) is valid.
   * The returned value is a reference to the coefficient of the internal
   * matrix, which may change if the matrix is modified.
   */
  inline const FT &get_mu_exp(int i, int j, long &expo);
  inline const FT &get_mu_exp(int i, int j);

  /**
   * Returns f = (b_i, b*_j) / ||b*_j||^2.
   *
   * Returns reference to `f`.
   */
  inline FT &get_mu(FT &f, int i, int j);

  /**
   * Return maximum bstar_i for all i
   */
  ZT get_max_gram();

  /**
   * Return maximum bstar_i for all i
   */
  FT get_max_bstar();

  /**
   * Returns f = r(i, j) and expo such that (b_i, b*_j) = f * 2^expo.
   * If enable_row_expo=false, expo is always 0.
   * If enable_row_expo=true, expo = row_expo[i] + row_expo[j]
   * If is assumed that r(i, j) is valid.
   * The returned value is a reference to the coefficient of the internal
   * matrix, which may change if the matrix is modified
   */
  inline const FT &get_r_exp(int i, int j, long &expo);
  inline const FT &get_r_exp(int i, int j);

  /**
   * Returns f = (b_i, b*_j).
   *
   * Returns reference to `f`.
   */
  inline FT &get_r(FT &f, int i, int j);

  /**
   * Returns expo such that mu(i, j) &lt; 2^expo for all j &lt; n_columns.
   * It is assumed that mu(i, j) is valid for all j &lt; n_columns.
   */
  long get_max_mu_exp(int i, int n_columns);

  /**
   * Updates r(i, j) and mu(i, j) if needed for all j in [0, last_j].
   * All coefficients of r and mu above the i-th row in columns
   * [0, min(last_j, i - 1)] must be valid.
   * If i=n_known_rows, n_known_rows is increased by one.
   */
  bool update_gso_row(int i, int last_j);

  /**
   * Updates r(i, j) and mu(i, j) if needed for all j.
   * All coefficients of r and mu above the i-th row in columns
   * [0, min(last_j, i - 1)] must be valid.
   * If i=n_known_rows, n_known_rows is increased by one.
   */
  inline bool update_gso_row(int i);

  /**
   * Updates all GSO coefficients (mu and r).
   */
  virtual inline bool update_gso();

  /**
   * Allows row_addmul(_we) for all rows even if the GSO has never been computed.
   */
  inline void discover_all_rows();

  /**
   * Sets the value of r(i, j). During the execution of LLL, some coefficients
   * are computed by the algorithm. They are set directly to avoid double
   * computation.
   */
  void set_r(int i, int j, FT &f);

  /**
   * Row old_r becomes row new_r and intermediate rows are shifted.
   * If new_r < old_r, then old_r must be < n_known_rows.
   */
  virtual void move_row(int old_r, int new_r) = 0;

  /**
   * b[i] := b[i] + x * b[j].
   * After one or several calls to row_addmul, row_op_end must be called.
   * Special cases |x| &lt;= 1 and |x| &lt;= LONG_MAX are optimized.
   * x should be an integer.
   * If row_op_force_long=true, x is always converted to (2^expo * long) instead
   * of (2^expo * ZT), which is faster if ZT=mpz_t but might lead to a loss of
   * precision (in LLL, more Babai iterations are needed).
   */
  virtual inline void row_addmul(int i, int j, const FT &x);

  /**
   * b[i] := b[i] + x * 2^expo_add * b[j].
   * After one or several calls to row_addmul_we, row_op_end must be called.
   * Special cases |x| &lt;= 1 and |x| &lt;= LONG_MAX are optimized.
   * x should be an integer.
   * If row_op_force_long=true, x is always converted to (2^expo * long) instead
   * of (2^expo * ZT), which is faster if ZT=mpz_t but might lead to a loss of
   * precision (in LLL, more Babai iterations are needed).
   */

  virtual void row_addmul_we(int i, int j, const FT &x, long expo_add) = 0;

  // b[i] += b[j] / b[i] -= b[j] (i > j)
  virtual void row_add(int i, int j) = 0;
  virtual void row_sub(int i, int j) = 0;

  /**
   * Early reduction
   * Allowed when enable_int_gram=false,
   * n_known_cols large enough to compute all the g(i,j)
   */
  void lock_cols();
  void unlock_cols();

  /**
   * Adds a zero row to b (and to u if enableTranform=true). One or several
   * operations can be performed on this row with row_addmul(_we), then
   * row_op_end must be called.
   * Do not use if enable_inverse_transform=true.
   */
  inline void create_row();
  virtual void create_rows(int n_new_rows) = 0;

  /**
   * Removes the last row of b (and of u if enable_transform=true).
   * Do not use if enable_inverse_transform=true.
   */
  inline void remove_last_row();
  virtual void remove_last_rows(int n_removed_rows) = 0;

  /**
   * Executes transformation by creating extra rows,
   * Calculating new entries, swapping the new rows with previous ones,
   * And then removing the excess rows
   */
  void apply_transform(const Matrix<FT> &transform, int src_base, int target_base);

  void apply_transform(const Matrix<FT> &transform, int src_base)
  {
    apply_transform(transform, src_base, src_base);
  }

  /**
   * Dump mu matrix to parameter `mu`.

   * When a double pointer is provided the caller must ensure it can hold
   * block_size^2 entries. When a vector is provided new entries are pushed to
   * the end. In particular, existing entries are not overwritten or cleared.
   *
   * @note No row discovery or update is performed prior to dumping the matrix.
   */

  inline void dump_mu_d(double *mu, int offset = 0, int block_size = -1);
  inline void dump_mu_d(vector<double> mu, int offset = 0, int block_size = -1);

  /**
   * Dump r vector to parameter `r`.

   * When a double pointer is provided the caller must ensure it can hold
   * block_size entries. When a vector is provided new entries are pushed to the
   * end. In particular, existing entries are not overwritten or cleared.
   *
   * @note No row discovery or update is performed prior to dumping the matrix.
   */

  inline void dump_r_d(double *r, int offset = 0, int block_size = -1);
  inline void dump_r_d(vector<double> &r, int offset = 0, int block_size = -1);

  /**
     @brief Return slope of the curve fitted to the lengths of the vectors from
     `start_row` to `stop_row`.

     The slope gives an indication of the quality of the basis.

     @param start_row start row (inclusive)
     @param stop_row  stop row (exclusive)
     @return
  */

  double get_current_slope(int start_row, int stop_row);

  /**
     @brief Return (squared) root determinant of the basis.

     @param start_row start row (inclusive)
     @param end_row   stop row (exclusive)
  */

  FT get_root_det(int start_row, int end_row);

  /**
     @brief Return log of the (squared) determinant of the basis.

     @param start_row start row (inclusive)
     @param end_row   stop row (exclusive)
  */

  FT get_log_det(int start_row, int end_row);

  /**
     @brief Return slide potential of the basis

     @param start_row  start row (inclusive)
     @param end_row    stop row (exclusive)
     @param block_size block size
  */

  FT get_slide_potential(int start_row, int end_row, int block_size);

  /**
   * Prints mu,r and g matrix to ostream os.
   *
   */
  inline void print_mu_r_g(ostream &os);

  /** Exact computation of dot products (i.e. with type ZT instead of FT) */
  const bool enable_int_gram;

  /** Normalization of each row of b by a power of 2. */
  const bool enable_row_expo;

  /** Computation of the transform matrix. */
  const bool enable_transform;

  /**
   * Computation of the inverse transform matrix (transposed).
   * This works only if enable_transform=true.
   * This matrix has very large coefficients, computing it is slow.
   */
  const bool enable_inverse_transform;

  /** Computation uses Householder transformation. */
  const bool enable_householder;

  /**
   * Changes the behaviour of row_addmul(_we).
   * See the description of row_addmul.
   */
  const bool row_op_force_long;

protected:
  /** Allocates matrices and arrays whose size depends on d (all but tmp_col_expo).
    * When enable_int_gram=false, initializes bf.
    */
  virtual void size_increased() = 0;

  /* Increases known rows and invalidates the last
   * gram row (gf) when enable_int_gram is false.
   * When enable_int_gram is true, it computes
   * the new inner products for the last gram
   * row (g).
   */
  virtual void discover_row() = 0;

  // Marks mu(i, j) and r(i, j) as invalid for j >= new_valid_cols
  void invalidate_gso_row(int i, int new_valid_cols = 0);
  /*  Upates the i-th row of bf. It does not invalidate anything, so the caller
    * must take into account that it might change row_expo.
    */
  virtual void update_bf(int i) = 0;
  /* Marks g(i, j) for all j <= i (but NOT for j > i) */
  virtual void invalidate_gram_row(int i) = 0;

  // b[i] <- b[i] + x * b[j] (i > j)
  virtual void row_addmul_si(int i, int j, long x) = 0;
  // b[i] <- b[i] + (2^expo * x) * b[j] (i > j)

  // TODO check if these must be scratched here!
  // virtual void row_addmul_si_2exp(int i, int j, long x, long expo) = 0;
  // virtual void row_addmul_2exp(int i, int j, const ZT &x, long expo) = 0;

  void symmetrize_g();

public:
  // Made public for bkz.cpp (dbkz_postprocessing)
  // b[i] <-> b[j] (i < j)
  virtual void row_swap(int i, int j) = 0;

protected:
  inline ZT &sym_g(int i, int j);

  /* Floating-point representation of the basis. It is used when
     enable_int_gram=false. */
  Matrix<FT> bf;

  Matrix<ZT> &u;        // Transform
  Matrix<ZT> &u_inv_t;  // Transposed inverse transform

  // init_row_size[i] = (last non-zero column in the i-th row of b) + 1
  vector<int> init_row_size;

  // bf[i], g[i], gf[i], mu[i] and r[i] are invalid for i >= n_known_rows
  int n_known_rows;
  int n_source_rows;  // Known rows before the beginning of early reduction
  int n_known_cols;
  bool cols_locked;
  int alloc_dim;

  /**
   * Coefficients of the Gram-Schmidt orthogonalization
   * (lower triangular matrix).
   *
   * If enable_row_expo=false,
   * mu(i, j) = (b_i, b*_j) / ||b*_j||^2.
   * If enable_row_expo=true,
   * mu(i, j) = (b_i, b*_j) / ||b*_j||^2  /  2 ^ (row_expo[i] - row_expo[j]).
   *
   * mu(i, j) is valid if 0 &lt;= i &lt; n_known_rows (&lt;= d) and
   * 0 &lt;= j &lt; min(gso_valid_cols[i], i)
   */
  Matrix<FT> mu;

  /**
   * Coefficients of the Gram-Schmidt orthogonalization
   * (lower triangular matrix).
   *
   * If enable_row_expo=false,
   * r(i, j) = (b_i, b*_j).
   * If enable_row_expo=true,
   * r(i, j) = (b_i, b*_j)  /  2 ^ (row_expo[i] + row_expo[j]).
   *
   * r(i, j) is valid if 0 &lt;= i &lt; n_known_rows (&lt;= d) and
   * 0 &lt;= j &lt; gso_valid_cols[i] (&lt;= i + 1).
   */
  Matrix<FT> r;

public:
  /** Replaced the gram matrix by a pointer. In the gso-class
    * there is also a matrix g, and in the constructor gptr is
    * defined to be equal to &g. In the ggso-class a gram matrix
    * is given (arg_g), and gptr is defined as &arg_g.
    */

  /* Gram matrix (dot products of basis vectors, lower triangular matrix)
   g(i, j) is valid if 0 <= i < n_known_rows and j <= i */
  Matrix<ZT> *gptr;
  //  Matrix<ZT> g;

protected:
  Matrix<FT> gf;

  /* Number of valid columns of the i-th row of mu and r.
     Valid only for 0 <= i < n_known_rows */
  vector<int> gso_valid_cols;

  /* Used by update_gso_row (+ update_gso), get_max_mu_exp and row_addmul_we. */
  FT ftmp1, ftmp2;
  /* Used by row_add, row_sub, row_addmul_si_2exp, row_addmul_2exp and
     indirectly by row_addmul. */
  ZT ztmp1;
  /* Used by row_addmul. */
  ZT ztmp2;
  /* Used by update_bf. */
  vector<long> tmp_col_expo;

#ifdef DEBUG
  /* Used only in debug mode. */
  int row_op_first, row_op_last;
  bool in_row_op_range(int i) { return i >= row_op_first && i < row_op_last; }
#endif
};

template <class ZT, class FT> inline MatGSOInterface<ZT, FT>::~MatGSOInterface()
{
  // delete gptr;
}

template <class ZT, class FT> inline void MatGSOInterface<ZT, FT>::symmetrize_g()
{
  if (gptr == nullptr)
  {
    throw std::runtime_error("Error: gptr is equal to the nullpointer.");
  }
  Matrix<ZT> &gr = *gptr;
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      gr(i, j) = sym_g(i, j);
    }
  }
}

template <class ZT, class FT> inline void MatGSOInterface<ZT, FT>::print_mu_r_g(ostream &os)
{
  os << "mu = " << endl;
  mu.print(os);
  os << endl << "r = " << endl;
  r.print(os);
  os << endl;
  if (gptr != nullptr)
  {
    os << "g = " << endl;
    symmetrize_g();
    gptr->print(os);
    os << endl << endl;
  }
}

template <class ZT, class FT> inline ZT &MatGSOInterface<ZT, FT>::sym_g(int i, int j)
{
  if (gptr == nullptr)
  {
    throw std::runtime_error("Error: gptr is equal to the nullpointer.");
  }
  Matrix<ZT> &gr = *gptr;
  return (i >= j) ? gr(i, j) : gr(j, i);
}

template <class ZT, class FT>
inline const FT &MatGSOInterface<ZT, FT>::get_mu_exp(int i, int j, long &expo)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && j >= 0 && j < i && j < gso_valid_cols[i] &&
                    !in_row_op_range(i));
  if (enable_row_expo)
    expo = row_expo[i] - row_expo[j];
  else
    expo = 0;
  return mu(i, j);
}

template <class ZT, class FT> inline const FT &MatGSOInterface<ZT, FT>::get_mu_exp(int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && j >= 0 && j < i && j < gso_valid_cols[i] &&
                    !in_row_op_range(i));
  return mu(i, j);
}

template <class ZT, class FT> inline FT &MatGSOInterface<ZT, FT>::get_mu(FT &f, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && j >= 0 && j < i && j < gso_valid_cols[i] &&
                    !in_row_op_range(i));

  f = mu(i, j);
  if (enable_row_expo)
    f.mul_2si(f, row_expo[i] - row_expo[j]);
  return f;
}

template <class ZT, class FT>
inline const FT &MatGSOInterface<ZT, FT>::get_r_exp(int i, int j, long &expo)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && j >= 0 && j < gso_valid_cols[i] &&
                    !in_row_op_range(i));
  if (enable_row_expo)
    expo = row_expo[i] + row_expo[j];
  else
    expo = 0;
  return r(i, j);
}

template <class ZT, class FT> inline const FT &MatGSOInterface<ZT, FT>::get_r_exp(int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && j >= 0 && j < gso_valid_cols[i] &&
                    !in_row_op_range(i));
  return r(i, j);
}

template <class ZT, class FT> inline FT &MatGSOInterface<ZT, FT>::get_r(FT &f, int i, int j)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && j >= 0 && j < gso_valid_cols[i] &&
                    !in_row_op_range(i));
  f = r(i, j);
  if (enable_row_expo)
    f.mul_2si(f, row_expo[i] + row_expo[j]);
  return f;
}

template <class ZT, class FT> inline bool MatGSOInterface<ZT, FT>::update_gso_row(int i)
{
  return update_gso_row(i, i);
}

template <class ZT, class FT> inline void MatGSOInterface<ZT, FT>::set_r(int i, int j, FT &f)
{
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && gso_valid_cols[i] >= j && j >= 0 &&
                    j < n_source_rows);
  r(i, j) = f;
  if (gso_valid_cols[i] == j)
    gso_valid_cols[i]++;
}

template <class ZT, class FT>
inline void MatGSOInterface<ZT, FT>::row_addmul(int i, int j, const FT &x)
{
  row_addmul_we(i, j, x, 0);
}

template <class ZT, class FT> inline void MatGSOInterface<ZT, FT>::create_row() { create_rows(1); }

template <class ZT, class FT> inline void MatGSOInterface<ZT, FT>::remove_last_row()
{
  remove_last_rows(1);
}

template <class ZT, class FT> inline void MatGSOInterface<ZT, FT>::discover_all_rows()
{
  while (n_known_rows < d)
    discover_row();
}

template <class ZT, class FT> inline bool MatGSOInterface<ZT, FT>::update_gso()
{
  for (int i = 0; i < d; i++)
  {
    if (!update_gso_row(i))
      return false;
  }
  return true;
}

#ifdef DEBUG
template <class ZT, class FT> inline void MatGSOInterface<ZT, FT>::row_op_begin(int first, int last)
{
  FPLLL_DEBUG_CHECK(row_op_first == -1);
  row_op_first = first;
  row_op_last  = last;
}
#else
template <class ZT, class FT>
inline void MatGSOInterface<ZT, FT>::row_op_begin(int /*first*/, int /*last*/)
{
}
#endif

template <class ZT, class FT>
inline void MatGSOInterface<ZT, FT>::dump_mu_d(double *mu, int offset, int block_size)
{
  FT e;
  if (block_size <= 0)
  {
    block_size = get_rows_of_b();
  }

  for (int i = 0; i < block_size; ++i)
  {
    for (int j = 0; j < block_size; ++j)
    {
      get_mu(e, offset + i, offset + j);
      mu[i * block_size + j] = e.get_d();
    }
  }
}

template <class ZT, class FT>
inline void MatGSOInterface<ZT, FT>::dump_mu_d(vector<double> mu, int offset, int block_size)
{
  FT e;
  if (block_size <= 0)
  {
    block_size = get_rows_of_b();
  }

  mu.reserve(mu.size() + block_size * block_size);
  for (int i = 0; i < block_size; ++i)
  {
    for (int j = 0; j < block_size; ++j)
    {
      get_mu(e, offset + i, offset + j);
      mu.push_back(e.get_d());
    }
  }
}

template <class ZT, class FT>
inline void MatGSOInterface<ZT, FT>::dump_r_d(double *r, int offset, int block_size)
{
  FT e;
  if (block_size <= 0)
  {
    block_size = get_rows_of_b();
  }

  for (int i = 0; i < block_size; ++i)
  {
    get_r(e, offset + i, offset + i);
    r[i] = e.get_d();
  }
}

template <class ZT, class FT>
inline void MatGSOInterface<ZT, FT>::dump_r_d(vector<double> &r, int offset, int block_size)
{
  FT e;
  if (block_size <= 0)
  {
    block_size = get_rows_of_b();
  }

  r.reserve(r.size() + block_size * block_size);
  for (int i = 0; i < block_size; ++i)
  {
    get_r(e, offset + i, offset + i);
    r.push_back(e.get_d());
  }
}

FPLLL_END_NAMESPACE

#endif
