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

#ifndef FPLLL_MATRIX_CPP
#define FPLLL_MATRIX_CPP

#include "matrix.h"
#include "../defs.h"

FPLLL_BEGIN_NAMESPACE

template <class T> void Matrix<T>::resize(int rows, int cols)
{
  int old_size = matrix.size();
  if (old_size < rows)
  {
    vector<NumVect<T>> m2(max(old_size * 2, rows));
    for (int i = 0; i < old_size; i++)
    {
      matrix[i].swap(m2[i]);
    }
    matrix.swap(m2);
  }
  for (int i = r; i < rows; i++)
  {
    matrix[i].resize(cols);
  }
  if (cols != c)
  {
    for (int i = min(r, rows) - 1; i >= 0; i--)
    {
      matrix[i].resize(cols);
    }
  }
  r    = rows;
  c    = cols;
  cols = c;
}

template <class T> template <class U> void Matrix<T>::fill(U value)
{
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < c; j++)
    {
      matrix[i][j] = value;
    }
  }
}

template <class T> void Matrix<T>::rotate_gram_left(int first, int last, int n_valid_rows)
{
  FPLLL_DEBUG_CHECK(0 <= first && first <= last && last < n_valid_rows && n_valid_rows <= r);
  matrix[first][first].swap(matrix[first][last]);
  for (int i = first; i < last; i++)
  {
    matrix[i + 1][first].swap(matrix[first][i]);
  }
  for (int i = first; i < n_valid_rows; i++)
  {
    matrix[i].rotate_left(first, min(last, i));  // most expensive step
  }
  rotate_left(first, last);
}

template <class T> void Matrix<T>::rotate_gram_right(int first, int last, int n_valid_rows)
{
  FPLLL_DEBUG_CHECK(0 <= first && first <= last && last < n_valid_rows && n_valid_rows <= r);
  rotate_right(first, last);
  for (int i = first; i < n_valid_rows; i++)
  {
    matrix[i].rotate_right(first, min(last, i));  // most expensive step
  }
  for (int i = first; i < last; i++)
  {
    matrix[i + 1][first].swap(matrix[first][i]);
  }
  matrix[first][first].swap(matrix[first][last]);
}

template <class T> void Matrix<T>::transpose()
{
  extend_vect(matrix, c);
  for (int i = 0; i < c; i++)
  {
    matrix[i].extend(r);
  }
  for (int i = 0; i < min(r, c); i++)
  {
    for (int j = i + 1; j < max(r, c); j++)
    {
      matrix[i][j].swap(matrix[j][i]);
    }
    if (c > r)
      matrix[i].resize(r);
  }
  std::swap(r, c);
}

template <class T> T Matrix<T>::get_max()
{
  T m, a;
  m = 0;
  for (int i = 0; i < r; i++)
    for (int j = 0; j < c; j++)
    {
      a.abs(matrix[i][j]);
      m = max(m, a);
    }
  return m;
}

template <class T> long Matrix<T>::get_max_exp()
{
  long max_exp = 0;
  for (int i = 0; i < r; i++)
    for (int j = 0; j < c; j++)
      max_exp = max(max_exp, matrix[i][j].exponent());
  return max_exp;
}

template <class T> void Matrix<T>::print(ostream &os, int nrows, int ncols) const
{
  if (nrows < 0 || nrows > r)
    nrows = r;
  if (ncols < 0 || ncols > c)
    ncols = c;
  os << '[';
  for (int i = 0; i < nrows; i++)
  {
    if (i > 0)
      os << '\n';
    os << '[';
    for (int j = 0; j < ncols; j++)
    {
      if (j > 0)
        os << ' ';
      os << matrix[i][j];
    }
    if (print_mode == MAT_PRINT_REGULAR && ncols > 0)
      os << ' ';
    os << ']';
  }
  if (print_mode == MAT_PRINT_REGULAR && nrows > 0)
    os << '\n';
  os << ']';
}

template <class T> void Matrix<T>::read(istream &is)
{
  char ch;
  matrix.clear();
  if (!(is >> ch))
    return;
  if (ch != '[')
  {
    is.setstate(ios::failbit);
    return;
  }
  while (is >> ch && ch != ']')
  {
    is.putback(ch);
    matrix.resize(matrix.size() + 1);
    if (!(is >> matrix.back()))
    {
      matrix.pop_back();
      break;
    }
  }

  r = matrix.size();
  c = 0;
  for (int i = 0; i < r; i++)
  {
    c = max(c, matrix[i].size());
  }
  for (int i = 0; i < r; i++)
  {
    int old_c = matrix[i].size();
    if (old_c < c)
    {
      matrix[i].resize(c);
      for (int j = old_c; j < c; j++)
      {
        matrix[i][j] = 0;
      }
    }
  }
}

/* ZZ_mat */

template <class ZT> inline void ZZ_mat<ZT>::gen_intrel(int bits)
{
  if (c != r + 1)
  {
    FPLLL_ABORT("gen_intrel called on an ill-formed matrix");
    return;
  }
  int i, j;
  for (i = 0; i < r; i++)
  {
    matrix[i][0].randb(bits);
    for (j = 1; j <= i; j++)
    {
      matrix[i][j] = 0;
    }
    matrix[i][i + 1] = 1;
    for (j = i + 2; j < c; j++)
    {
      matrix[i][j] = 0;
    }
  }
}

template <class ZT> inline void ZZ_mat<ZT>::gen_simdioph(int bits, int bits2)
{
  if (c != r)
  {
    FPLLL_ABORT("gen_simdioph called on an ill-formed matrix");
    return;
  }
  int i, j;

  matrix[0][0] = 1;
  matrix[0][0].mul_2si(matrix[0][0], bits2);
  for (i = 1; i < r; i++)
    matrix[0][i].randb(bits);
  for (i = 1; i < r; i++)
  {
    for (j         = 1; j < i; j++)
      matrix[j][i] = 0;
    matrix[i][i]   = 1;
    matrix[i][i].mul_2si(matrix[i][i], bits);
    for (j         = i + 1; j < c; j++)
      matrix[j][i] = 0;
  }
}

template <class ZT> inline void ZZ_mat<ZT>::gen_uniform(int bits)
{
  if (c != r)
  {
    FPLLL_ABORT("gen_uniform called on an ill-formed matrix");
    return;
  }
  for (int i = 0; i < r; i++)
    for (int j = 0; j < c; j++)
      matrix[i][j].randb(bits);
}

template <class ZT> inline void ZZ_mat<ZT>::gen_ntrulike(int bits)
{
  // [A00 A01]
  // [A10 A11]

  int i, j, k;
  int d = r / 2;
  if (c != r || c != 2 * d)
  {
    FPLLL_ABORT("gen_ntrulike called on an ill-formed matrix");
    return;
  }
  // clang-format off
  Z_NR<ZT> *h = new Z_NR<ZT>[d];
  // clang-format on
  Z_NR<ZT> q;

  q.randb(bits);
  if (q.sgn() == 0)
    q  = 1;
  h[0] = 0;
  for (i = 1; i < d; i++)
  {
    h[i].randm(q);
    h[0].sub(h[0], h[i]);
    if (h[0].sgn() < 0)
      h[0].add(h[0], q);
    // set h0 such that h(1) = 0 mod q.
  }

  // I in A00
  for (i = 0; i < d; i++)
  {
    for (j         = 0; j < i; j++)
      matrix[i][j] = 0;
    matrix[i][i]   = 1;
    for (j         = i + 1; j < d; j++)
      matrix[i][j] = 0;
  }

  // 0 in A10
  for (i = d; i < r; i++)
  {
    for (j         = 0; j < d; j++)
      matrix[i][j] = 0;
  }
  // qI in A11
  for (i = d; i < r; i++)
  {
    for (j         = d; j < i; j++)
      matrix[i][j] = 0;
    matrix[i][i]   = q;
    for (j         = i + 1; j < c; j++)
      matrix[i][j] = 0;
  }
  // H in A01
  for (i = 0; i < d; i++)
    for (j = d; j < c; j++)
    {
      k = j - d - i;
      while (k < 0)
      {
        k += d;
      }
      matrix[i][j] = h[k];
    }

  delete[] h;
}

template <class ZT> inline void ZZ_mat<ZT>::gen_ntrulike_withq(int q)
{
  // Same as above, except q is specified by the user rather than chosen
  // randomly with a prescribed bit-length.

  int i, j, k;
  int d = r / 2;
  if (c != r || c != 2 * d)
  {
    FPLLL_ABORT("gen_ntrulike called on an ill-formed matrix");
    return;
  }
  // clang-format off
  Z_NR<ZT> *h = new Z_NR<ZT>[d];
  // clang-format on
  Z_NR<ZT> q2;

  q2   = q;
  h[0] = 0;
  for (i = 1; i < d; i++)
  {
    h[i].randm(q2);
    h[0].sub(h[0], h[i]);
    if (h[0].sgn() < 0)
      h[0].add(h[0], q2);
    // set h0 such that h(1) = 0 mod q.
  }

  // I in A00
  for (i = 0; i < d; i++)
  {
    for (j         = 0; j < i; j++)
      matrix[i][j] = 0;
    matrix[i][i]   = 1;
    for (j         = i + 1; j < d; j++)
      matrix[i][j] = 0;
  }

  // 0 in A10
  for (i = d; i < r; i++)
  {
    for (j         = 0; j < d; j++)
      matrix[i][j] = 0;
  }
  // qI in A11
  for (i = d; i < r; i++)
  {
    for (j         = d; j < i; j++)
      matrix[i][j] = 0;
    matrix[i][i]   = q2;
    for (j         = i + 1; j < c; j++)
      matrix[i][j] = 0;
  }
  // H in A01
  for (i = 0; i < d; i++)
    for (j = d; j < c; j++)
    {
      k = j - d - i;
      while (k < 0)
      {
        k += d;
      }
      matrix[i][j] = h[k];
    }

  delete[] h;
}

template <class ZT> inline void ZZ_mat<ZT>::gen_ntrulike2(int bits)
{

  int i, j, k;
  int d = r / 2;
  if (c != r || c != 2 * d)
  {
    FPLLL_ABORT("gen_ntrulike2 called on an ill-formed matrix");
    return;
  }
  // clang-format off
  Z_NR<ZT> *h = new Z_NR<ZT>[d];
  Z_NR<ZT> q;
  // clang-format on

  q.randb(bits);
  h[0] = 0;
  for (i = 1; i < d; i++)
  {
    h[i].randm(q);
    h[0].sub(h[0], h[i]);
    if (h[0].sgn() < 0)
      h[0].add(h[0], q);
    // set h0 such that h(1) = 0 mod q.
  }

  for (i = 0; i < d; i++)
  {
    for (j         = 0; j < c; j++)
      matrix[i][j] = 0;
  }

  for (i         = 0; i < d; i++)
    matrix[i][i] = q;
  for (i = d; i < r; i++)
    for (j         = d; j < c; j++)
      matrix[i][j] = 0;
  for (i         = d; i < c; i++)
    matrix[i][i] = 1;

  for (i = d; i < r; i++)
  {
    for (j = 0; j < d; j++)
    {
      k = i - d - j;
      while (k < 0)
      {
        k += d;
      }
      matrix[i][j] = h[k];
    }
  }

  delete[] h;
}

template <class ZT> inline void ZZ_mat<ZT>::gen_ntrulike2_withq(int q)
{

  int i, j, k;
  int d = r / 2;
  if (c != r || c != 2 * d)
  {
    FPLLL_ABORT("gen_ntrulike2 called on an ill-formed matrix");
    return;
  }
  // clang-format off
  Z_NR<ZT> *h = new Z_NR<ZT>[d];
  Z_NR<ZT> q2;
  // clang-format on

  q2   = q;
  h[0] = 0;
  for (i = 1; i < d; i++)
  {
    h[i].randm(q2);
    h[0].sub(h[0], h[i]);
    if (h[0].sgn() < 0)
      h[0].add(h[0], q2);
    // set h0 such that h(1) = 0 mod q.
  }

  for (i = 0; i < d; i++)
  {
    for (j         = 0; j < c; j++)
      matrix[i][j] = 0;
  }

  for (i         = 0; i < d; i++)
    matrix[i][i] = q2;
  for (i = d; i < r; i++)
    for (j         = d; j < c; j++)
      matrix[i][j] = 0;
  for (i         = d; i < c; i++)
    matrix[i][i] = 1;

  for (i = d; i < r; i++)
  {
    for (j = 0; j < d; j++)
    {
      k = i - d - j;
      while (k < 0)
      {
        k += d;
      }
      matrix[i][j] = h[k];
    }
  }

  delete[] h;
}

template <class ZT> inline void ZZ_mat<ZT>::gen_qary(int k, Z_NR<ZT> &q)
{
  int i, j;
  int d = r;
  if (c != r || k > r)
  {
    FPLLL_ABORT("gen_qary called on an ill-formed matrix");
    return;
  }

  for (i = 0; i < d - k; i++)
    for (j         = 0; j < d - k; j++)
      matrix[i][j] = 0;

  for (i         = 0; i < d - k; i++)
    matrix[i][i] = 1;

  for (i = 0; i < d - k; i++)
    for (j = d - k; j < d; j++)
      matrix[i][j].randm(q);

  for (i = d - k; i < d; i++)
    for (j         = 0; j < d - k; j++)
      matrix[i][j] = 0;

  for (i         = d - k; i < d; i++)
    matrix[i][i] = q;
}

template <class ZT> inline void ZZ_mat<ZT>::gen_trg(double alpha)
{
  int i, j, bits;
  Z_NR<ZT> ztmp, ztmp2, zone, sign;

  ztmp2 = 0;
  zone  = 1;

  int d = r;
  if (c != r)
  {
    FPLLL_ABORT("gen_trg called on an ill-formed matrix");
    return;
  }

  for (i = 0; i < d; i++)
  {
    bits = (int)pow((double)(2 * d - i), alpha);
    ztmp = 1;
    ztmp.mul_2si(ztmp, bits);
    ztmp.sub(ztmp, zone);
    matrix[i][i].randm(ztmp);
    matrix[i][i].add_ui(matrix[i][i], 2);
    ztmp.div_2si(matrix[i][i], 1);
    for (j = i + 1; j < d; j++)
    {
      matrix[j][i].randm(ztmp);
      sign.randb(1);
      if (sign == 1)
        matrix[j][i].sub(ztmp2, matrix[j][i]);
      matrix[i][j] = 0;
    }
  }
}

template <class ZT> inline void ZZ_mat<ZT>::gen_trg2(FP_NR<mpfr_t> *w)
{
  int i, j;
  Z_NR<ZT> ztmp, ztmp2;

  int d = r;
  if (c != r)
  {
    FPLLL_ABORT("gen_trg2 called on an ill-formed matrix");
    return;
  }

  for (i = 0; i < d; i++)
  {
    matrix[i][i].set_f(w[i]);
    ztmp.div_2si(matrix[i][i], 1);
    ztmp2 = 1;
    ztmp.add(ztmp, ztmp2);
    for (j = i + 1; j < d; j++)
    {
      ztmp2 = 0;
      matrix[j][i].randm(ztmp);
      if (rand() % 2 == 1)
        matrix[j][i].sub(ztmp2, matrix[j][i]);
      matrix[i][j] = 0;
    }
  }
}

FPLLL_END_NAMESPACE

#endif
