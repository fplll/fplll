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

#ifndef FPLLL_GSO_H
#define FPLLL_GSO_H

#include "nr/matrix.h"

FPLLL_BEGIN_NAMESPACE

enum MatGSOFlags {
  GSO_DEFAULT = 0,
  GSO_INT_GRAM = 1,
  GSO_ROW_EXPO = 2,
  GSO_OP_FORCE_LONG = 4
};

/**
 * MatGSO provides an interface for performing elementary operations on a basis
 * and computing its Gram matrix and its Gram-Schmidt orthogonalization.
 * The Gram-Schmidt coefficients are computed on demand. The object keeps track
 * of which coefficients are valid after each row operation.
 */
template<class ZT, class FT>
class MatGSO {
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
   * @param uInvT
   *   Inverse transform (should be empty, which disables the computation, or
   *   initialized with identity matrix). It works only if u is not empty.
   * @param enableIntGram
   *   If true, coefficients of the Gram matrix are computed with exact integer
   *   arithmetic (type ZT). Otherwise, they are computed in floating-point
   *   (type FT). Note that when exact arithmetic is used, all coefficients of
   *   the first nKnownRows are continuously updated, whereas in floating-point,
   *   they are computed only on-demand. This option cannot be enabled if
   *   enableRowExpo=true.
   * @param enableRowExpo
   *   If true, each row of b is normalized by a power of 2 before doing
   *   conversion to floating-point, which hopefully avoids some overflows.
   *   This option cannot be enabled if enableIntGram=true and works only
   *   with FT=double and FT=long double. It is useless and MUST NOT be used
   *   for FT=dpe or FT=mpfr_t.
   * @param rowOpForceLong
   *   Affects the behaviour of row_addmul(_we).
   *   See the documentation of row_addmul.
   */
  MatGSO(Matrix<ZT>& b, Matrix<ZT>& u, Matrix<ZT>& uInvT, int flags);

  /**
   * Number of rows of b (dimension of the lattice).
   * Can be changed with createRow or removeLastRow.
   */
  int d;

  /**
   * Basis of the lattice
   */
  Matrix<ZT>& b;

  /**
   * When enableRowExpo=true, rowExpo[i] is the smallest non-negative integer
   * such that b(i, j) &lt;= 2^rowExpo[i] for all j. Otherwise this array is empty.
   */
  vector<long> rowExpo;

  /**
   * Must be called before a sequence of row_addmul(_we).
   */
  inline void rowOpBegin(int first, int last);

  /**
   * Must be called after a sequence of row_addmul(_we). This invalidates the
   * i-th line of the GSO.
   */
  void rowOpEnd(int first, int last);

  /**
   * Returns Gram matrix coefficients (0 &lt;= i &lt; nKnownRows and
   * 0 &lt;= j &lt;= i).
   * If enableRowExpo=false, returns the dot product (b[i], b[j]).
   * If enableRowExpo=true, returns
   * (b[i], b[j]) / 2 ^ (rowExpo[i] + rowExpo[j]).
   */
  inline void getGram(FT& f, int i, int j);

  /**
   * Returns the mu matrix
   * Coefficients of the Gram Schmidt Orthogonalization
   * (lower triangular matrix)
   * mu(i, j) = r(i, j) / ||b*_j||^2.
   */
  const Matrix<FT>& getMuMatrix() {
    return mu;
  }
  
    /**
   * Returns the r matrix
   * Coefficients of the Gram Schmidt Orthogonalization
   * (lower triangular matrix)
   */
  const Matrix<FT>& getRMatrix() {
    return r;
  }

  /**
   * Returns f = mu(i, j) and expo such that
   * f * 2^expo = (b_i, b*_j) / ||b*_j||^2.
   * If enableRowExpo=false, expo is always 0.
   * If enableRowExpo=true, expo = rowExpo[i] - rowExpo[j]
   * It is assumed that mu(i, j) is valid.
   * The returned value is a reference to the coefficient of the internal
   * matrix, which may change if the matrix is modified.
   */
  inline const FT& getMuExp(int i, int j, long& expo);
  inline const FT& getMuExp(int i, int j);

  /**
   * Returns f = (b_i, b*_j) / ||b*_j||^2.
   */
  inline void getMu(FT& f, int i, int j);

  /**
   * Returns f = r(i, j) and expo such that (b_i, b*_j) = f * 2^expo.
   * If enableRowExpo=false, expo is always 0.
   * If enableRowExpo=true, expo = rowExpo[i] + rowExpo[j]
   * If is assumed that r(i, j) is valid.
   * The returned value is a reference to the coefficient of the internal
   * matrix, which may change if the matrix is modified
   */
  inline const FT& getRExp(int i, int j, long& expo);
  inline const FT& getRExp(int i, int j);

  /**
   * Returns f = (b_i, b*_j).
   */
  inline void getR(FT& f, int i, int j);

  /** 
   * Returns expo such that mu(i, j) &lt; 2^expo for all j &lt; nColumns.
   * It is assumed that mu(i, j) is valid for all j &lt; nColumns.
   */
  long getMaxMuExp(int i, int nColumns);

  /**
   * Updates r(i, j) and mu(i, j) if needed for all j in [0, lastJ].
   * All coefficients of r and mu above the i-th row in columns
   * [0, min(lastJ, i - 1)] must be valid.
   * If i=nKnownRows, nKnownRows is increased by one.
   */
  bool updateGSORow(int i, int lastJ);

  /**
   * Updates r(i, j) and mu(i, j) if needed for all j.
   * All coefficients of r and mu above the i-th row in columns
   * [0, min(lastJ, i - 1)] must be valid.
   * If i=nKnownRows, nKnownRows is increased by one.
   */
  inline bool updateGSORow(int i);

  /**
   * Updates all GSO coefficients (mu and r).
   */
  inline bool updateGSO();

  /**
   * Allows row_addmul(_we) for all rows even if the GSO has never been computed.
   */
  inline void discoverAllRows();

  /**
   * Sets the value of r(i, j). During the execution of LLL, some coefficients
   * are computed by the algorithm. They are set directly to avoid double
   * computation.
   */
  void setR(int i, int j, FT& f);

  /**
   * Row oldR becomes row newR and intermediate rows are shifted.
   * If newR < oldR, then oldR must be < nKnownRows.
   */
  void moveRow(int oldR, int newR);

  /**
   * b[i] := b[i] + x * b[j].
   * After one or several calls to row_addmul, rowOpEnd must be called.
   * Special cases |x| &lt;= 1 and |x| &lt;= LONG_MAX are optimized.
   * x should be an integer.
   * If rowOpForceLong=true, x is always converted to (2^expo * long) instead
   * of (2^expo * ZT), which is faster if ZT=mpz_t but might lead to a loss of
   * precision (in LLL, more Babai iterations are needed).
   */
  inline void row_addmul(int i, int j, const FT& x);

  /**
   * b[i] := b[i] + x * 2^expoAdd * b[j].
   * After one or several calls to row_addmul_we, rowOpEnd must be called.
   * Special cases |x| &lt;= 1 and |x| &lt;= LONG_MAX are optimized.
   * x should be an integer.
   * If rowOpForceLong=true, x is always converted to (2^expo * long) instead
   * of (2^expo * ZT), which is faster if ZT=mpz_t but might lead to a loss of
   * precision (in LLL, more Babai iterations are needed).
   */
  void row_addmul_we(int i, int j, const FT& x, long expoAdd);

  /** 
   * Early reduction
   * Allowed when enableIntGram=false,
   * nKnownCols large enough to compute all the g(i,j)
   */
  void lockCols();
  void unlockCols();

  /**
   * Adds a zero row to b (and to u if enableTranform=true). One or several
   * operations can be performed on this row with row_addmul(_we), then
   * rowOpEnd must be called.
   * Do not use if enableInvTransform=true.
   */
  inline void createRow();
  inline void createRows(int nNewRows);

  /**
   * Removes the last row of b (and of u if enableTransform=true).
   * Do not use if enableInvTransform=true.
   */
  inline void removeLastRow();
  inline void removeLastRows(int nRemovedRows);

  /**
   * Executes transformation by creating extra rows,
   * Calculating new entries, swapping the new rows with previous ones,
   * And then removing the excess rows
   */
  void applyTransform(const Matrix<FT>& transform, int srcBase, int targetBase);

  void applyTransform(const Matrix<FT>& transform, int srcBase) {
    applyTransform(transform, srcBase, srcBase);
  }

  /**
   * Dump mu matrix to parameter `mu`.

   * When a double pointer is provided the caller must ensure it can hold
   * blocksize^2 entries. When a vector is provided new entries are pushed to
   * the end. In particular, existing entries are not overwritten or cleared.
   *
   * @note No row discovery or update is performed prior to dumping the matrix.
   */

  inline void dumpMu_d(double* mu, int offset=0, int blocksize=-1);
  inline void dumpMu_d(vector<double> mu, int offset=0, int blocksize=-1);

  /**
   * Dump r vector to parameter `r`.

   * When a double pointer is provided the caller must ensure it can hold
   * blocksize entries. When a vector is provided new entries are pushed to the
   * end. In particular, existing entries are not overwritten or cleared.
   *
   * @note No row discovery or update is performed prior to dumping the matrix.
   */

  inline void dumpR_d(double* r, int offset=0, int blocksize=-1);
  inline void dumpR_d(vector<double> r, int offset=0, int blocksize=-1);

  /** Exact computation of dot products (i.e. with type ZT instead of FT) */
  const bool enableIntGram;

  /** Normalization of each row of b by a power of 2. */
  const bool enableRowExpo;

  /** Computation of the transform matrix. */
  const bool enableTransform;

  /**
   * Computation of the inverse transform matrix (transposed).
   * This works only if enableTransform=true.
   * This matrix has very large coefficients, computing it is slow.
   */
  const bool enableInvTransform;

  /**
   * Changes the behaviour of row_addmul(_we).
   * See the description of row_addmul.
   */
  const bool rowOpForceLong;

private:
  /* Allocates matrices and arrays whose size depends on d (all but tmpColExpo).
     When enableIntGram=false, initializes bf. */
  void sizeIncreased();

  void discoverRow();

  // Marks mu(i, j) and r(i, j) as invalid for j >= newValidCols
  inline void invalidateGSORow(int i, int newValidCols = 0);
  /* Upates the i-th row of bf. It does not invalidate anything, so the caller
     must take into account that it might change rowExpo. */
  void updateBF(int i);
  /* Marks g(i, j) for all j <= i (but NOT for j > i) */
  void invalidateGramRow(int i);

  // b[i] += b[j] / b[i] -= b[j] (i > j)
  void row_add(int i, int j);
  void row_sub(int i, int j);
  // b[i] <- b[i] + x * b[j] (i > j)
  void row_addmul_si(int i, int j, long x);
  // b[i] <- b[i] + (2^expo * x) * b[j] (i > j)
  void row_addmul_si_2exp(int i, int j, long x, long expo);
  void row_addmul_2exp(int i, int j, const ZT& x, long expo);
  // b[i] <-> b[j] (i < j)
  void rowSwap(int i, int j);

  inline ZT& symG(int i, int j) {
    return (i >= j) ? g(i, j) : g(j, i);
  }

  /* Floating-point representation of the basis. It is used when
     enableIntGram=true. */
  Matrix<FT> bf;

  Matrix<ZT>& u;        // Transform
  Matrix<ZT>& uInvT;    // Transposed inverse transform

  // initRowSize[i] = (last non-zero column in the i-th row of b) + 1
  vector<int> initRowSize;

  // bf[i], g[i], gf[i], mu[i] and r[i] are invalid for i >= nKnownRows
  int nKnownRows;
  int nSourceRows; // Known rows before the beginning of early reduction
  int nKnownCols;
  bool colsLocked;
  int allocDim;

  /**
   * Coefficients of the Gram-Schmidt orthogonalization
   * (lower triangular matrix).
   *
   * If enableRowExpo=false,
   * mu(i, j) = (b_i, b*_j) / ||b*_j||^2.
   * If enableRowExpo=true,
   * mu(i, j) = (b_i, b*_j) / ||b*_j||^2  /  2 ^ (rowExpo[i] - rowExpo[j]).
   *
   * mu(i, j) is valid if 0 &lt;= i &lt; nKnownRows (&lt;= d) and
   * 0 &lt;= j &lt; min(gsoValidCols[i], i)
   */
  Matrix<FT> mu;

  /**
   * Coefficients of the Gram-Schmidt orthogonalization
   * (lower triangular matrix).
   *
   * If enableRowExpo=false,
   * r(i, j) = (b_i, b*_j).
   * If enableRowExpo=true,
   * r(i, j) = (b_i, b*_j)  /  2 ^ (rowExpo[i] + rowExpo[j]).
   *
   * r(i, j) is valid if 0 &lt;= i &lt; nKnownRows (&lt;= d) and
   * 0 &lt;= j &lt; gsoValidCols[i] (&lt;= i + 1).
   */
  Matrix<FT> r;

  /* Gram matrix (dot products of basis vectors, lower triangular matrix)
     g(i, j) is valid if 0 <= i < nKnownRows and j <= i */
  Matrix<ZT> g;
  Matrix<FT> gf;

  /* Number of valid columns of the i-th row of mu and r.
     Valid only for 0 <= i < nKnownRows */
  vector<int> gsoValidCols;

  /* Used by updateGSORow (+ updateGSO), getMaxMuExp and row_addmul_we. */
  FT ftmp1, ftmp2;
  /* Used by row_add, row_sub, row_addmul_si_2exp, row_addmul_2exp and
     indirectly by row_addmul. */
  ZT ztmp1;
  /* Used by row_addmul. */
  ZT ztmp2;
  /* Used by updateBF. */
  vector<long> tmpColExpo;

#ifdef DEBUG
  /* Used only in debug mode. */
  int rowOpFirst, rowOpLast;
  bool inRowOpRange(int i) {
    return i >= rowOpFirst && i < rowOpLast;
  }
#endif
};

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::getGram(FT& f, int i, int j) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows && j >= 0 && j <= i
                    && j < nSourceRows && !inRowOpRange(i));
  if (enableIntGram)
    f.set_z(g(i, j));
  else {
    if (gf(i, j).is_nan()) {
      dotProduct(gf(i, j), bf[i], bf[j], nKnownCols);
    }
    f = gf(i, j);
  }
}

template<class ZT, class FT>
inline const FT& MatGSO<ZT, FT>::getMuExp(int i, int j, long& expo) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows && j >= 0 && j < i
                    && j < gsoValidCols[i] && !inRowOpRange(i));
  if (enableRowExpo)
    expo = rowExpo[i] - rowExpo[j];
  else
    expo = 0;
  return mu(i, j);
}

template<class ZT, class FT>
inline const FT& MatGSO<ZT, FT>::getMuExp(int i, int j) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows && j >= 0 && j < i
                    && j < gsoValidCols[i] && !inRowOpRange(i));
  return mu(i, j);
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::getMu(FT& f, int i, int j) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows && j >= 0 && j < i
                    && j < gsoValidCols[i] && !inRowOpRange(i));
  f = mu(i, j);
  if (enableRowExpo)
    f.mul_2si(f, rowExpo[i] - rowExpo[j]);
}

template<class ZT, class FT>
inline const FT& MatGSO<ZT, FT>::getRExp(int i, int j, long& expo) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows && j >= 0
                    && j < gsoValidCols[i] && !inRowOpRange(i));
  if (enableRowExpo)
    expo = rowExpo[i] + rowExpo[j];
  else
    expo = 0;
  return r(i, j);
}

template<class ZT, class FT>
inline const FT& MatGSO<ZT, FT>::getRExp(int i, int j) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows && j >= 0
                    && j < gsoValidCols[i] && !inRowOpRange(i));
  return r(i, j);
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::getR(FT& f, int i, int j) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows && j >= 0
                    && j < gsoValidCols[i] && !inRowOpRange(i));
  f = r(i, j);
  if (enableRowExpo)
    f.mul_2si(f, rowExpo[i] + rowExpo[j]);
}

template<class ZT, class FT>
inline bool MatGSO<ZT, FT>::updateGSORow(int i) {
  return updateGSORow(i, i);
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::setR(int i, int j, FT& f) {
  FPLLL_DEBUG_CHECK(i >= 0 && i < nKnownRows && gsoValidCols[i] >= j
                    && j >= 0 && j < nSourceRows);
  r(i, j) = f;
  if (gsoValidCols[i] == j) gsoValidCols[i]++;
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::row_addmul(int i, int j, const FT& x) {
  row_addmul_we(i, j, x, 0);
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::createRow() {
  createRows(1);
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::removeLastRow() {
  removeLastRows(1);
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::createRows(int nNewRows) {
  FPLLL_DEBUG_CHECK(!colsLocked);
  int oldD = d;
  d += nNewRows;
  b.setRows(d);
  for (int i = oldD; i < d; i++) {
    for (int j = 0; j < b.getCols(); j++) {
      b[i][j] = 0;
    }
  }
  if (enableTransform) {
    u.setRows(d);
    for (int i = oldD; i < d; i++)
      for (int j = 0; j < u.getCols(); j++)
        u[i][j] = 0;
  }
  sizeIncreased();
  if (nKnownRows == oldD) discoverAllRows();
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::removeLastRows(int nRemovedRows) {
  FPLLL_DEBUG_CHECK(!colsLocked && d >= nRemovedRows);
  d -= nRemovedRows;
  nKnownRows = min(nKnownRows, d);
  nSourceRows = nKnownRows;
  b.setRows(d);
  if (enableTransform)
    u.setRows(d);
}


template<class ZT, class FT>
inline void MatGSO<ZT, FT>::discoverAllRows() {
  while (nKnownRows < d)
    discoverRow();
}

template<class ZT, class FT>
inline bool MatGSO<ZT, FT>::updateGSO() {
  for (int i = 0; i < d; i++) {
    if (!updateGSORow(i))
      return false;
  }
  return true;
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::rowOpBegin(int first, int last) {
#ifdef DEBUG
  FPLLL_DEBUG_CHECK(rowOpFirst == -1);
  rowOpFirst = first;
  rowOpLast = last;
#endif
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::dumpMu_d(double* mu, int offset, int blocksize){
  FT e;
  if (blocksize <= 0) {
    blocksize = b.getRows();
  }

  for (int i = 0; i < blocksize; ++i) {
    for (int j = 0; j < blocksize; ++j) {
      getMu(e,offset+i,offset+j);
      mu[i*blocksize+j] = e.get_d();
    }
  }
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::dumpMu_d(vector<double> mu, int offset, int blocksize){
  FT e;
  if (blocksize <= 0) {
    blocksize = b.getRows();
  }

  r.reserve(r.size() + blocksize*blocksize);
  for (int i = 0; i < blocksize; ++i) {
    for (int j = 0; j < blocksize; ++j) {
      getMu(e,offset+i,offset+j);
      mu.push_back(e.get_d());
    }
  }
}


template<class ZT, class FT>
inline void MatGSO<ZT, FT>::dumpR_d(double* r, int offset, int blocksize){
  FT e;
  if (blocksize <= 0) {
    blocksize = b.getRows();
  }

  for (int i = 0; i < blocksize; ++i) {
    getR(e,offset+i, offset+i);
    r[i] = e.get_d();
  }
}

template<class ZT, class FT>
inline void MatGSO<ZT, FT>::dumpR_d(vector<double> r, int offset, int blocksize){
  FT e;
  if (blocksize <= 0) {
    blocksize = b.getRows();
  }

  r.reserve(r.size() + blocksize*blocksize);
  for (int i = 0; i < blocksize; ++i) {
    getR(e,offset+i,offset+i);
    r.push_back(e.get_d());
  }
}


FPLLL_END_NAMESPACE



#endif
