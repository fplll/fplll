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

#ifndef FPLLL_MATRIX_H
#define FPLLL_MATRIX_H

#include "numvect.h"

FPLLL_BEGIN_NAMESPACE

enum MatPrintMode {
  MAT_PRINT_COMPACT = 0,
  MAT_PRINT_REGULAR = 1
};

template<class T>
class Matrix;

/** MatrixRow stores a reference to a row of a Matrix. It supports a subset
    of operations available on vectors. */
template<class T>
class MatrixRow {
public:
  /** Returns a reference to the i-th element of this row. */
  T& operator[](int i)              {return row[i];}
  /** Returns a const reference to the i-th element of this row on constant
      objects. */
  const T& operator[](int i) const  {return row[i];}
  /** Returns the number of columns. */
  int size() const                  {return row.size();}
  /** Prints this object on stream os. */
  void print(ostream& os) const     {os << row;}

  bool is_zero(int from = 0) const {
    return row.is_zero(from);
  }
  int sizeNZ() const {
    return row.sizeNZ();
  }
  void fill(long value) {
    row.fill(value);
  }
  void add(const MatrixRow<T>& v) {
    row.add(v.row);
  }
  void add(const MatrixRow<T>& v, int n) {
    row.add(v.row, n);
  }
  void sub(const MatrixRow<T>& v) {
    row.sub(v.row);
  }
  void sub(const MatrixRow<T>& v, int n) {
    row.sub(v.row, n);
  }
  void addmul_2exp(const MatrixRow<T>& v, const T& x, long expo, T& tmp) {
    row.addmul_2exp(v.row, x, expo, tmp);
  }
  void addmul_2exp(const MatrixRow<T>& v, const T& x, long expo, int n, T& tmp) {
    row.addmul_2exp(v.row, x, expo, tmp);
  }
  void addmul_si(const MatrixRow<T>& v, long x) {
    row.addmul_si(v.row, x);
  }
  void addmul_si(const MatrixRow<T>& v, long x, int n) {
    row.addmul_si(v.row, x, n);
  }
  void addmul_si_2exp(const MatrixRow<T>& v, long x, long expo, T& tmp) {
    row.addmul_si_2exp(v.row, x, expo, tmp);
  }
  void addmul_si_2exp(const MatrixRow<T>& v, long x, long expo, int n, T& tmp) {
    row.addmul_si_2exp(v.row, x, expo, n, tmp);
  }

  friend class Matrix<T>;
private:
  MatrixRow(const NumVect<T>& row) : row(const_cast<NumVect<T>&>(row)) {}
  NumVect<T>& row;
};

template<class T>
void dotProduct(T& result, const MatrixRow<T>& v1,
                       const MatrixRow<T>& v2, int n) {
  FPLLL_DEBUG_CHECK(n > 0 && n <= v1.size() && v1.size() == v2.size()
                    && (v1.is_zero(n) || v2.is_zero(n)));
  result.mul(v1[0], v2[0]);
  for (int i = 1; i < n; i++) {
    result.addmul(v1[i], v2[i]);
  }
}

template<class T>
inline void dotProduct(T& result, const MatrixRow<T>& v1, const MatrixRow<T>& v2) {
  dotProduct(result, v1, v2, v1.size());
}

/** Prints a MatrixRow on stream os. */
template<class T>
ostream& operator<<(ostream& os, const MatrixRow<T>& row) {
  row.print(os);
  return os;
}

/** Matrix is a two-dimensional container. Read and write operations on single
    elements are in constant time. The amortized complexity of resizing the
    matrix is proportional to the number of added/removed elements. All indices
    are 0-based. */
template<class T>
class Matrix {
public:
  /** Creates an empty matrix (0 x 0). */
  Matrix() : r(0), c(0) {}
  /** Creates a matrix of dimensions rows x cols. All elements are
      initialized with the default constructor of T. */
  Matrix(int rows, int cols) : r(0), c(0) {
    resize(rows, cols);
  }

  /** Sets number of rows and the number of columns to 0. */
  void clear() {
    r = c = 0;
    matrix.clear();
  }
  /** Returns true if the matrix has 0 rows, false otherwise. */
  bool empty() const {
    return r == 0;
  }
  /** Sets the dimensions of this matrix, preserving as much as possible of the
      content. The value of new elements is undefined. */
  void resize(int rows, int cols);
  /** Sets the number of rows. Content is not erased except for deleted rows.
      The value of new elements is undefined. */
  void setRows(int rows) {
    resize(rows, c);
  }
  /** Sets the number of columns. Content is not erased except for deleted
      columns. The value of new elements is undefined. */
  void setCols(int cols) {
    resize(r, cols);
  }
  /** Fills this matrix with a given value. */
  template<class U>
  void fill(U value);
  /** Efficiently swaps two matrices. */
  void swap(Matrix<T>& m) {
    matrix.swap(m.matrix);
    std::swap(r, m.r);
    std::swap(c, m.c);
  }

  /** Returns the number of rows */
  int getRows() const {
    return r;
  }
  /** Returns the number of columns */
  int getCols() const {
    return c;
  }
  /** Returns a reference to the element (i, j). */
  T& operator()(int i, int j) {
    FPLLL_DEBUG_CHECK(i >= 0 && i < r && j >= 0 && j < c);
    return matrix[i][j];
  }
  /** Returns a constant reference to the element (i, j) on constant objects. */
  const T& operator()(int i, int j) const {
    FPLLL_DEBUG_CHECK(i >= 0 && i < r && j >= 0 && j < c);
    return matrix[i][j];
  }
  /** Returns a MatrixRow object pointing to the i-th row of this matrix. */
  MatrixRow<T> operator[](int i) {
    FPLLL_DEBUG_CHECK(i >= 0 && i < r);
    return MatrixRow<T>(matrix[i]);
  }
  /** Returns a MatrixRow object pointing to the i-th row of this matrix
      on constant objects. */
  const MatrixRow<T> operator[](int i) const {
    FPLLL_DEBUG_CHECK(i >= 0 && i < r);
    return MatrixRow<T>(matrix[i]);
  }
  /** Rows swap. */
  void swapRows(int r1, int r2) {
    matrix[r1].swap(matrix[r2]);
  }
  /** Rows permutation.
      (m[first],...,m[last]) becomes (m[first+1],...,m[last],m[first]) */
  void rotateLeft(int first, int last) {
    rotateLeftBySwap(matrix, first, last);
  }
  /** Rows permutation.
      (m[first],...,m[last]) becomes (m[last],m[first],...,m[last-1]) */
  void rotateRight(int first, int last) {
    rotateRightBySwap(matrix, first, last);
  }
  /** Rows permutation.
      (m[first],...,m[middle-1],m[middle],m[last]) becomes
      (m[middle],...,m[last],m[first],...,m[middle-1]) */
  void rotate(int first, int middle, int last) {
    rotateBySwap(matrix, first, middle, last);
  }
  /** Transformation needed to update the lower triangular Gram matrix when
     rotateLeft(first, last) is done on the basis of the lattice. */
  void rotateGramLeft(int first, int last, int nValidRows);
  /** Transformation needed to update the lower triangular Gram matrix when
      rotateRight(first, last) is done on the basis of the lattice. */
  void rotateGramRight(int first, int last, int nValidRows);
  /** Transpose. */
  void transpose();
  long getMaxExp();
  /** Prints this matrix. No end-of-line is printed after the last line.
      @param os      output stream
      @param nRows   maximum number of rows to display (optional)
      @param nCols   maximum number of columns to display (optional) */
  void print(ostream& os, int nRows = -1, int nCols = -1) const;
  /** Reads this matrix from a stream. */
  void read(istream& is);

#ifdef FPLLL_V3_COMPAT
  // Old interface (do not use)
  int GetNumCols() const          {return c;}
  int GetNumRows() const          {return r;}
  void SetNumCols(int cols)       {resize(r, cols);}
  void SetNumRows(int rows)       {resize(rows, c);}
  T& Get(int i, int j)            {return matrix[i][j];}
  void Set(int i, int j, T& s)    {matrix[i][j] = s;}
  T* GetVec(int i)                {return &matrix[i][0];}
  const T* GetVecC(int i) const   {return &matrix[i][0];}
  void print(int d, int n);
  void print()                    {print(r, c);}
  int read();
#endif

  static int setPrintMode(int newPrintMode) {
    int oldMode = printMode;
    printMode = newPrintMode;
    return oldMode;
  }

protected:
  int r, c;
  vector<NumVect<T> > matrix;

  static int printMode;
};

template<class T>
int Matrix<T>::printMode = MAT_PRINT_COMPACT;

template<class T>
ostream& operator<<(ostream& os, const Matrix<T>& m) {
  m.print(os);
  return os;
}

template<class T>
istream& operator>>(istream& is, Matrix<T>& m) {
  m.read(is);
  return is;
}

/** ZZ_mat is a matrix of integers. */
template <class ZT>
class ZZ_mat : public Matrix<Z_NR<ZT> > {
public:
  typedef Z_NR<ZT> T;
  using Matrix<T>::r;
  using Matrix<T>::c;
  using Matrix<T>::matrix;
  using Matrix<T>::resize;
  using Matrix<T>::getCols;
  using Matrix<T>::getRows;

  /** Creates an empty matrix (0 x 0). */
  ZZ_mat() : Matrix<T>() {}
  /** Creates a matrix of dimensions rows x cols. All elements are
      initialized with the default constructor of Z_NR&lt;T&gt;. */
  ZZ_mat(int rows, int cols) : Matrix<T>(rows, cols) {}

  // generators
  void gen_zero(int d, int n) {
    resize(d, n);
    for (int i = 0; i < d; i++)
      matrix[i].fill(0);
  }
  void gen_identity(int d) {
    gen_zero(d, d);
    for (int i = 0; i < d; i++)
      matrix[i][i] = 1;
  }
  void gen_intrel(int bits);
  void gen_simdioph(int bits,int bits2);
  void gen_uniform(int bits);


  /** Construct a matrix `[[I,H],[0,qI]]` where `H` is constructed from rotations of a vector ``h``.

     @note The constructed matrix will not come with a guarantee of unusually short vectors.
  **/

  void gen_ntrulike(int bits);
  void gen_ntrulike_withq(int q);

  /** Construct a matrix ``[[qI,0],[H,I]]`` where ``H`` is constructed from rotations of a vector ``h``.

    @note The constructed matrix will not come with a guarantee of unusually short vectors.
  */

  void gen_ntrulike2(int bits);
  void gen_ntrulike2_withq(int q);

  /** Construct a matrix ``[[qI,0],[H,I]]`` where ``H`` is uniform mod q, of dimensions k x (n-k).
  */

  void gen_qary(int k, int bits);
  void gen_qary_withq(int k, int q);

  /** Construct a lower triangular matrices with specified diagonal coefficients
and random sub-diagonal coefficients.
  */

  void gen_trg(double alpha);
  void gen_trg2(FP_NR<mpfr_t> *w);

#ifdef FPLLL_V3_COMPAT
  int getShift();
#endif
};

/** FP_mat is a matrix of floating-point numbers. */
template <class FT>
class FP_mat : public Matrix<FP_NR<FT> > {
public:
  typedef FP_NR<FT> T;
  /** Creates an empty matrix (0 x 0). */
  FP_mat() : Matrix<T>() {}
  /** Creates a matrix of dimensions rows x cols. All elements are
      initialized with the default constructor of FP_NR&lt;T&gt;. */
  FP_mat(int rows, int cols) : Matrix<T>(rows, cols) {}
};

typedef ZZ_mat<IntegerT> IntMatrix;
typedef FP_mat<FloatT> FloatMatrix;

FPLLL_END_NAMESPACE

#include "matrix.cpp"

#endif
