/* Copyright (C) 2011 Xavier Pujol.

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

#ifndef FPLLL_NUMVECT_H
#define FPLLL_NUMVECT_H

#include "fplll/nr/nr.h"
#include <vector>

FPLLL_BEGIN_NAMESPACE

/** Extends the size of the given vector. */
template <class T> void extend_vect(vector<T> &v, int size)
{
  if (static_cast<int>(v.size()) < size)
  {
    v.resize(size);
  }
}

/** Reverses a vector by consecutive swaps. */
template <class T> void reverse_by_swap(vector<T> &v, int first, int last)
{
  for (; first < last; first++, last--)
    v[first].swap(v[last]);
}

/** Rotates a vector by consecutive swaps. */
template <class T> void rotate_by_swap(vector<T> &v, int first, int middle, int last)
{
  // Algorithm from STL code
  reverse_by_swap(v, first, middle - 1);
  reverse_by_swap(v, middle, last);
  for (; first < middle && middle <= last; first++, last--)
  {
    v[first].swap(v[last]);
  }
  if (first == middle)
    reverse_by_swap(v, middle, last);
  else
    reverse_by_swap(v, first, middle - 1);
}

/** Rotates a vector left-wise by consecutive swaps. */
template <class T> void rotate_left_by_swap(vector<T> &v, int first, int last)
{
  FPLLL_DEBUG_CHECK(0 <= first && first <= last && last < static_cast<int>(v.size()));
  for (int i = first; i < last; i++)
  {
    v[i].swap(v[i + 1]);
  }
}

/** Rotates a vector right-wise by consecutive swaps. */
template <class T> void rotate_right_by_swap(vector<T> &v, int first, int last)
{
  FPLLL_DEBUG_CHECK(0 <= first && first <= last && last < static_cast<int>(v.size()));
  for (int i = last - 1; i >= first; i--)
  {
    v[i].swap(v[i + 1]);
  }
}

/** Print a vector on stream os. */
template <class T> ostream &operator<<(ostream &os, const vector<T> &v)
{
  os << "[";
  int n = v.size();
  for (int i = 0; i < n; i++)
  {
    if (i > 0)
      os << " ";
    os << v[i];
  }
  os << "]";
  return os;
}

/** Reads a vector from stream is. */
template <class T> istream &operator>>(istream &is, vector<T> &v)
{
  char c;
  v.clear();
  if (!(is >> c))
    return is;
  if (c != '[')
  {
    is.setstate(ios::failbit);
    return is;
  }
  while (is >> c && c != ']')
  {
    is.putback(c);
    v.resize(v.size() + 1);
    if (!(is >> v.back()))
    {
      v.pop_back();
      return is;
    }
  }
  return is;
}

/** Generate a zero vector. */
template <class T> void gen_zero_vect(vector<T> &v, int n)
{
  v.resize(n);
  fill(v.begin(), v.end(), 0);
}

template <class T> class NumVect;

template <class T> ostream &operator<<(ostream &os, const NumVect<T> &v);

template <class T> istream &operator>>(istream &is, NumVect<T> &v);

template <class T> class NumVect
{
public:
  typedef typename vector<T>::iterator iterator;
  /** Creates an empty NumVect (0). */
  NumVect() {}
  /** Initializes NumVect with the elements of the given NumVect. */
  NumVect(const NumVect &v) : data(v.data) {}
  /** Initializes NumVect of specific size, The initial content is
      undefined. */
  NumVect(int size) : data(size) {}
  NumVect(int size, const T &t) : data(size, t) {}
  /** Sets the NumVect to the elements of the given NumVect. */
  void operator=(const NumVect &v)
  {
    if (this != &v)
      data = v.data;
  }
  /** Swaps the data between the NumVect and the given NumVect. */
  void swap(NumVect &v) { data.swap(v.data); }
  /** Returns an iterator to the beginning of NumVect.*/
  const iterator begin() { return data.begin(); }
  /** Returns an iterator to the end of NumVect. */
  iterator end() { return data.end(); }
  /** Returns the number of elements in NumVect. */
  int size() const { return data.size(); }
  /** Checks whether NumVect is empty. */
  bool empty() const { return data.empty(); }
  /** Sets the size of NumVect. */
  void resize(int size) { data.resize(size); }
  /** Sets the size of NumVect and all its elemnts to t.*/
  void resize(int size, const T &t) { data.resize(size, t); }
  /** Sets the size of NumVect and all its elements to zero. */
  void gen_zero(int size)
  {
    data.resize(size);
    fill(0);
  }
  /** Inserts an element in the back of NumVect. */
  void push_back(const T &t) { data.push_back(t); }
  /** Removes the back element of NumVect. */
  void pop_back() { data.pop_back(); }
  /** Returns a reference to the front element of NumVect. */
  T &front() { return data.front(); }
  /** Returns a const reference to the front element of NumVect,
      on constant object. */
  const T &front() const { return data.front(); }
  /** Returns a reference to the back element of NumVect. */
  T &back() { return data.back(); }
  /** Returns a const reference to the back element of NumVect,
      on constant object. */
  const T &back() const { return data.back(); }
  /** Extends the size of NumVect, only if it is needed. */
  void extend(int maxSize)
  {
    if (size() < maxSize)
      data.resize(maxSize);
  }
  void clear() { data.clear(); }
  /** Returns a reference to the i-th element of NumVect. */
  T &operator[](int i) { return data[i]; }
  /** Returns a const reference to the i-th element of NumVect on constant
      object. */
  const T &operator[](int i) const { return data[i]; }
  /** Addition of two NumVector objects, till index n. */
  void add(const NumVect<T> &v, int n);
  /** Addition of two NumVector objects. */
  void add(const NumVect<T> &v) { add(v, size()); }
  /** Subtraction of two NumVector objects, till index n. */
  void sub(const NumVect<T> &v, int n);
  /** Subtraction of two NumVector objects. */
  void sub(const NumVect<T> &v) { sub(v, size()); }
  /** Multiplication of NumVector and a number c, from index b till index n. */
  void mul(const NumVect<T> &v, int b, int n, T c);
  /** Multiplication of NumVector and a number c, till index n. */
  void mul(const NumVect<T> &v, int n, T c);
  /** Multiplication of NumVector and a number c. */
  void mul(const NumVect<T> &v, T c) { mul(v, size(), c); }
  /** Division of NumVector and a number c, from index b till index n. */
  void div(const NumVect<T> &v, int b, int n, T c);
  /** Division of NumVector and a number c, till index n. */
  void div(const NumVect<T> &v, int n, T c);
  /** Division of NumVector and a number c. */
  void div(const NumVect<T> &v, T c) { div(v, size(), c); }
  /** Incremeanting each coefficient of NumVector by its product with
      number c, from beg to index n - 1. */
  void addmul(const NumVect<T> &v, T x, int beg, int n);
  /** Incremeanting each coefficient of NumVector by its product with
      number c, till index n. */
  void addmul(const NumVect<T> &v, T x, int n);
  /** Incremeanting each coefficient of NumVector by its product with
      number c. */
  void addmul(const NumVect<T> &v, T x) { addmul(v, x, size()); }
  void addmul_2exp(const NumVect<T> &v, const T &x, long expo, T &tmp)
  {
    addmul_2exp(v, x, expo, size(), tmp);
  }
  void addmul_2exp(const NumVect<T> &v, const T &x, long expo, int n, T &tmp);
  void addmul_si(const NumVect<T> &v, long x) { addmul_si(v, x, size()); }
  void addmul_si(const NumVect<T> &v, long x, int n);
  void addmul_si_2exp(const NumVect<T> &v, long x, long expo, T &tmp)
  {
    addmul_si_2exp(v, x, expo, size(), tmp);
  }
  void addmul_si_2exp(const NumVect<T> &v, long x, long expo, int n, T &tmp);

  /** (v[first],...,v[last]) becomes (v[first+1],...,v[last],v[first]) */
  void rotate_left(int first, int last) { rotate_left_by_swap(data, first, last); }

  /** (v[first],...,v[last]) becomes (v[last],v[first],...,v[last-1]) */
  void rotate_right(int first, int last) { rotate_right_by_swap(data, first, last); }

  /** Returns expo >= 0 such that all elements are < 2^expo. */
  long get_max_exponent();

  /** Fills NumVect with the value given. */
  void fill(long value);

  /** Checks if NumVect has zero elements from fromCol. */
  bool is_zero(int fromCol = 0) const;

  /** Returns last non-zero index of NumVector. */
  int size_nz() const;

  friend ostream &operator<<<T>(ostream &os, const NumVect<T> &v);
  friend istream &operator>><T>(istream &is, NumVect<T> &v);

private:
  vector<T> data;
};

template <class T> void NumVect<T>::add(const NumVect<T> &v, int n)
{
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
    data[i].add(data[i], v[i]);
}

template <class T> void NumVect<T>::sub(const NumVect<T> &v, int n)
{
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
    data[i].sub(data[i], v[i]);
}

template <class T> void NumVect<T>::mul(const NumVect<T> &v, int b, int n, T c)
{
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= b; i--)
    data[i].mul(v[i], c);
}

template <class T> void NumVect<T>::mul(const NumVect<T> &v, int n, T c) { mul(v, 0, n, c); }

template <class T> void NumVect<T>::div(const NumVect<T> &v, int b, int n, T c)
{
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= b; i--)
    data[i].div(v[i], c);
}

template <class T> void NumVect<T>::div(const NumVect<T> &v, int n, T c) { div(v, 0, n, c); }

template <class T> void NumVect<T>::addmul(const NumVect<T> &v, T x, int beg, int n)
{
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= beg; i--)
    data[i].addmul(v[i], x);
}

template <class T> void NumVect<T>::addmul(const NumVect<T> &v, T x, int n)
{
  this->addmul(v, x, 0, n);
}

template <class T>
void NumVect<T>::addmul_2exp(const NumVect<T> &v, const T &x, long expo, int n, T &tmp)
{
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
  {
    tmp.mul(v[i], x);
    tmp.mul_2si(tmp, expo);
    data[i].add(data[i], tmp);
  }
}

template <class T> void NumVect<T>::addmul_si(const NumVect<T> &v, long x, int n)
{
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
    data[i].addmul_si(v[i], x);
}

template <class T>
void NumVect<T>::addmul_si_2exp(const NumVect<T> &v, long x, long expo, int n, T &tmp)
{
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
  {
    tmp.mul_si(v[i], x);
    tmp.mul_2si(tmp, expo);
    data[i].add(data[i], tmp);
  }
}

template <class T> long NumVect<T>::get_max_exponent()
{
  long max_expo = 0;
  for (int i = 0; i < size(); i++)
  {
    max_expo = max(max_expo, data[i].exponent());
  }
  return max_expo;
}

template <class T> void NumVect<T>::fill(long value)
{
  for (int i = 0; i < size(); i++)
  {
    data[i] = value;
  }
}

template <class T> bool NumVect<T>::is_zero(int fromCol) const
{
  for (int i = fromCol; i < size(); i++)
  {
    if (!data[i].is_zero())
      return false;
  }
  return true;
}

template <class T> int NumVect<T>::size_nz() const
{
  int i;
  for (i = data.size(); i > 0; i--)
  {
    if (data[i - 1] != 0)
      break;
  }
  return i;
}

/** Compute the truncated dot product between tow Numvect using coefficients [beg, n).
 * Constraint: n > beg.
 */
template <class T>
inline void dot_product(T &result, const NumVect<T> &v1, const NumVect<T> &v2, int beg, int n)
{
  FPLLL_DEBUG_CHECK(beg >= 0 && n > beg && n <= v1.size() && n <= v2.size());
  //(v1.is_zero(n) || v2.is_zero(n))); tested previously
  result.mul(v1[beg], v2[beg]);
  for (int i = beg + 1; i < n; i++)
  {
    result.addmul(v1[i], v2[i]);
  }
}

template <class T>
inline void dot_product(T &result, const NumVect<T> &v1, const NumVect<T> &v2, int n)
{
  FPLLL_DEBUG_CHECK(n <= v1.size() && v1.size() == v2.size() && (v1.is_zero(n) || v2.is_zero(n)));
  dot_product(result, v1, v2, 0, n);
}

template <class T> inline void dot_product(T &result, const NumVect<T> &v1, const NumVect<T> &v2)
{
  dot_product(result, v1, v2, v1.size());
}

template <class T> inline void squared_norm(T &result, const NumVect<T> &v)
{
  dot_product(result, v, v);
}

/** Prints a NumVect on stream os. */
template <class T> ostream &operator<<(ostream &os, const NumVect<T> &v) { return os << v.data; }

/** Reads a NumVect from stream is. */
template <class T> istream &operator>>(istream &is, NumVect<T> &v) { return is >> v.data; }

FPLLL_END_NAMESPACE

#endif
