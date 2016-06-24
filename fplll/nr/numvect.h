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

#include "nr.h"
#include <vector>

FPLLL_BEGIN_NAMESPACE

template<class T>
void extendVect(vector<T>& v, int size) {
  if (static_cast<int>(v.size()) < size) {
    v.resize(size);
  }
}

template<class T>
void reverseBySwap(vector<T>& v, int first, int last) {
  for (; first < last; first++, last--)
    v[first].swap(v[last]);
}

template<class T>
void rotateBySwap(vector<T>& v, int first, int middle, int last) {
  // Algorithm from STL code
  reverseBySwap(v, first, middle - 1);
  reverseBySwap(v, middle, last);
  for (; first < middle && middle <= last; first++, last--) {
    v[first].swap(v[last]);
  }
  if (first == middle)
    reverseBySwap(v, middle, last);
  else
    reverseBySwap(v, first, middle - 1);
}

template<class T>
void rotateLeftBySwap(vector<T>& v, int first, int last) {
  FPLLL_DEBUG_CHECK(0 <= first && first <= last
          && last < static_cast<int>(v.size()));
  for (int i = first; i < last; i++) {
    v[i].swap(v[i + 1]);
  }
}

template<class T>
void rotateRightBySwap(vector<T>& v, int first, int last) {
  FPLLL_DEBUG_CHECK(0 <= first && first <= last
          && last < static_cast<int>(v.size()));
  for (int i = last - 1; i >= first; i--) {
    v[i].swap(v[i + 1]);
  }
}

/** Prints a vector on stream os. */
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v) {
  os << "[";
  int n = v.size();
  for (int i = 0; i < n; i++) {
    if (i > 0) os << " ";
    os << v[i];
  }
  os << "]";
  return os;
}

/** Reads a vector from stream is. */
template<class T>
istream& operator>>(istream& is, vector<T>& v) {
  char c;
  v.clear();
  if (!(is >> c)) return is;
  if (c != '[') {
    is.setstate(ios::failbit);
    return is;
  }
  while (is >> c && c != ']') {
    is.putback(c);
    v.resize(v.size() + 1);
    if (!(is >> v.back())) {
      v.pop_back();
      return is;
    }
  }
  return is;
}

template<class T>
void genZeroVect(vector<T>& v, int n) {
  v.resize(n);
  fill(v.begin(), v.end(), 0);
}



template<class T>
class NumVect;

template<class T>
ostream& operator<<(ostream& os, const NumVect<T>& v);

template<class T>
istream& operator>>(istream& is, NumVect<T>& v);



template<class T>
class NumVect {
public:
  typedef typename vector<T>::iterator iterator;

  NumVect() {}
  NumVect(const NumVect& v) : data(v.data) {}
  NumVect(int size) : data(size) {} // Initial content is undefined
  NumVect(int size, const T& t) : data(size, t) {}

  void operator=(const NumVect& v) {
    if (this != &v)
      data = v.data;
  }
  void swap(NumVect& v) {
    data.swap(v.data);
  }

  const iterator begin() {
    return data.begin();
  }
  iterator end() {
    return data.end();
  }
  int size() const {
    return data.size();
  }
  bool empty() const {
    return data.empty();
  }
  void resize(int size) {
    data.resize(size);
  }
  void resize(int size, const T& t) {
    data.resize(size, t);
  }
  void gen_zero(int size) {
    data.resize(size);
    fill(0);
  }
  void push_back(const T& t) {
    data.push_back(t);
  }
  void pop_back() {
    data.pop_back();
  }
  T& front() {
    return data.front();
  }
  const T& front() const {
    return data.front();
  }
  T& back() {
    return data.back();
  }
  const T& back() const {
    return data.back();
  }
  void extend(int maxSize) {
    if (size() < maxSize)
      data.resize(maxSize);
  }
  void clear() {
    data.clear();
  }
  T& operator[](int i) {
    return data[i];
  }
  const T& operator[](int i) const {
    return data[i];
  }

  void add(const NumVect<T>& v, int n);
  void add(const NumVect<T>& v) {add(v, size());}
  void sub(const NumVect<T>& v, int n);
  void sub(const NumVect<T>& v) {sub(v, size());}
  void mul(const NumVect<T>& v, int n, T c);
  void mul(const NumVect<T>& v, T c) {mul(v, size(), c);}
  void addmul(const NumVect<T>& v, T x, int n);
  void addmul(const NumVect<T>& v, T x) {addmul(v, x, size());}
  void addmul_2exp(const NumVect<T>& v, const T& x, long expo, T& tmp) {
    addmul_2exp(v, x, expo, size(), tmp);
  }
  void addmul_2exp(const NumVect<T>& v, const T& x, long expo, int n, T& tmp);
  void addmul_si(const NumVect<T>& v, long x) {addmul_si(v, x, size());}
  void addmul_si(const NumVect<T>& v, long x, int n);
  void addmul_si_2exp(const NumVect<T>& v, long x,
                             long expo, T& tmp) {
    addmul_si_2exp(v, x, expo, size(), tmp);
  }
  void addmul_si_2exp(const NumVect<T>& v, long x,
                             long expo, int n, T& tmp);

  /** (v[first],...,v[last]) becomes (v[first+1],...,v[last],v[first]) */
  void rotateLeft(int first, int last) {
    rotateLeftBySwap(data, first, last);
  }

  /** (v[first],...,v[last]) becomes (v[last],v[first],...,v[last-1]) */
  void rotateRight(int first, int last) {
    rotateRightBySwap(data, first, last);
  }

  /** Returns expo >= 0 such that all elements are < 2^expo. */
  long getMaxExponent();

  void fill(long value);

  bool is_zero(int fromCol = 0) const;

  int sizeNZ() const;

  friend ostream& operator<< <T>(ostream& os, const NumVect<T>& v);
  friend istream& operator>> <T>(istream& is, NumVect<T>& v);

private:
  vector<T> data;
};

template<class T>
void NumVect<T>::add(const NumVect<T>& v, int n) {
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
    data[i].add(data[i], v[i]);
}

template<class T>
void NumVect<T>::sub(const NumVect<T>& v, int n) {
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
    data[i].sub(data[i], v[i]);
}

template<class T>
void NumVect<T>::mul(const NumVect<T>& v, int n, T c) {
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
    data[i].mul(v[i], c);
}

template<class T>
void NumVect<T>::addmul(const NumVect<T>& v, T x, int n) {
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
    data[i].addmul(v[i], x);
}

template<class T>
void NumVect<T>::addmul_2exp(const NumVect<T>& v, const T& x,
                        long expo, int n, T& tmp) {
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--) {
    tmp.mul(v[i], x);
    tmp.mul_2si(tmp, expo);
    data[i].add(data[i], tmp);
  }
}

template<class T>
void NumVect<T>::addmul_si(const NumVect<T>& v, long x, int n) {
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--)
    data[i].addmul_si(v[i], x);
}

template<class T>
void NumVect<T>::addmul_si_2exp(const NumVect<T>& v, long x,
                            long expo, int n, T& tmp) {
  FPLLL_DEBUG_CHECK(n <= size() && size() == v.size() && v.is_zero(n));
  for (int i = n - 1; i >= 0; i--) {
    tmp.mul_si(v[i], x);
    tmp.mul_2si(tmp, expo);
    data[i].add(data[i], tmp);
  }
}

template<class T>
long NumVect<T>::getMaxExponent() {
  long maxExpo = 0;
  for (int i = 0; i < size(); i++) {
    maxExpo = max(maxExpo, data[i].exponent());
  }
  return maxExpo;
}

template<class T>
void NumVect<T>::fill(long value) {
  for (int i = 0; i < size(); i++) {
    data[i] = value;
  }
}

template<class T>
bool NumVect<T>::is_zero(int fromCol) const {
  for (int i = fromCol; i < size(); i++) {
    if (!data[i].is_zero()) return false;
  }
  return true;
}

template<class T>
int NumVect<T>::sizeNZ() const {
  int i;
  for (i = data.size(); i > 0; i--) {
    if (data[i - 1] != 0) break;
  }
  return i;
}




template<class T>
void dotProduct(T& result, const NumVect<T>& v1,
                       const NumVect<T>& v2, int n) {
  FPLLL_DEBUG_CHECK(n > 0 && n <= v1.size() && v1.size() == v2.size()
                    && (v1.is_zero(n) || v2.is_zero(n)));
  result.mul(v1[0], v2[0]);
  for (int i = 1; i < n; i++) {
    result.addmul(v1[i], v2[i]);
  }
}

template<class T>
inline void dotProduct(T& result, const NumVect<T>& v1, const NumVect<T>& v2) {
  dotProduct(result, v1, v2, v1.size());
}

template<class T>
inline void squaredNorm(T& result, const NumVect<T>& v) {
  dotProduct(result, v, v);
}

/** Prints a NumVect on stream os. */
template<class T>
ostream& operator<<(ostream& os, const NumVect<T>& v) {
  return os << v.data;
}

/** Reads a NumVect from stream is. */
template<class T>
istream& operator>>(istream& is, NumVect<T>& v) {
  return is >> v.data;
}

typedef vector<Integer> IntVect;
typedef vector<Float> FloatVect;

FPLLL_END_NAMESPACE

#endif
