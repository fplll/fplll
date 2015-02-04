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

#ifndef FPLLL_NR_CPP
#define FPLLL_NR_CPP

FPLLL_BEGIN_NAMESPACE

/**********************
 * int specialization *
 **********************/

#ifdef FPLLL_WITH_ZLONG

inline long computeLongExponent(long data) {
  unsigned long y = static_cast<unsigned long>(abs(data));
  long e;
  for (e = 0; y; e++, y >>= 1) {}
  return e;
}

template<>
inline Z_NR<long>::Z_NR() {}
template<>
inline Z_NR<long>::Z_NR(const Z_NR<long>& z) : data(z.data) {}
template<>
inline Z_NR<long>::~Z_NR() {}

template<>
inline double Z_NR<long>::get_d() const {
  return static_cast<double>(data);
}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<>
inline long double Z_NR<long>::get_ld() const {
  return static_cast<long double>(data);
}
#endif
template<>
inline long Z_NR<long>::get_si() const {
  return data;
}
template<>
inline void Z_NR<long>::set_str(const char* s) {
  data = atol(s);
}
template<>
inline void Z_NR<long>::operator=(const Z_NR<long>& a) {
  data = a.data;
}
template<>
inline void Z_NR<long>::operator=(long a) {
  data = static_cast<long>(a);
}

template<>
inline int Z_NR<long>::sgn() const {
  if (data > 0) return 1;
  if (data == 0) return 0;
  return -1;
}
template<>
inline int Z_NR<long>::cmp(const Z_NR<long>& m) const {
  if (data > m.data) return 1;
  if (data == m.data) return 0;
  return -1;
}
template<>
inline bool Z_NR<long>::operator<(const Z_NR<long>& a) const {
  return data < a.data;
}
template<>
inline bool Z_NR<long>::operator<(long a) const {
  return data < a;
}
template<>
inline bool Z_NR<long>::operator>(const Z_NR<long>& a) const {
  return data > a.data;
}
template<>
inline bool Z_NR<long>::operator>(long a) const {
  return data > a;
}
template<>
inline bool Z_NR<long>::operator<=(const Z_NR<long>& a) const {
  return data <= a.data;
}
template<>
inline bool Z_NR<long>::operator<=(long a) const {
  return data <= a;
}
template<>
inline bool Z_NR<long>::operator>=(const Z_NR<long>& a) const {
  return data >= a.data;
}
template<>
inline bool Z_NR<long>::operator>=(long a) const {
  return data >= a;
}
template<>
inline bool Z_NR<long>::operator==(const Z_NR<long>& a) const {
  return data == a.data;
}
template<>
inline bool Z_NR<long>::operator==(long a) const {
  return data == a;
}
template<>
inline bool Z_NR<long>::operator!=(const Z_NR<long>& a) const {
  return data != a.data;
}
template<>
inline bool Z_NR<long>::operator!=(long a) const {
  return data != a;
}

template<>
inline void Z_NR<long>::add(const Z_NR<long>& a, const Z_NR<long>& b) {
  data = a.data + b.data;
}
template<>
inline void Z_NR<long>::add_ui(const Z_NR<long>& a, unsigned int b) {
  data = a.data + b;
}
template<>
inline void Z_NR<long>::sub(const Z_NR<long>& a, const Z_NR<long>& b) {
  data = a.data - b.data;
}
template<>
inline void Z_NR<long>::neg(const Z_NR<long>& a) {
  data = -a.data;
}
template<>
inline void Z_NR<long>::mul(const Z_NR<long>& a, const Z_NR<long>& b) {
  data = a.data * b.data;
}
template<>
inline void Z_NR<long>::mul_si(const Z_NR<long>& a, long b) {
  data = a.data * b;
}
template<>
inline void Z_NR<long>::mul_ui(const Z_NR<long>& a, unsigned long b) {
  data = a.data * b;
}
template<>
inline void Z_NR<long>::mul_2si(const Z_NR<long>& a, long b) {
  //NOTE: if b >= size_in_bits(a), the result is undefined
  if (b >= 0)
    data = a.data << b;
  else
    data = a.data >> -b;
}
template<>
inline void Z_NR<long>::div_2si(const Z_NR<long>& a, long b) {
  //NOTE: if b >= size_in_bits(a), the result is undefined
  if (b >= 0)
    data = a.data >> b;
  else
    data = a.data << -b;
}

template<>
inline void Z_NR<long>::addmul(const Z_NR<long>& a, const Z_NR<long>& b) {
  data += a.data * b.data;
}
template<>
inline void Z_NR<long>::addmul_ui(const Z_NR<long>& a, unsigned long b) {
  data += a.data * b;
}
template<>
inline void Z_NR<long>::addmul_si(const Z_NR<long>& a, long b) {
  data += a.data * b;
}
template<>
inline void Z_NR<long>::submul(const Z_NR<long>& a, const Z_NR<long>& b) {
  data -= a.data * b.data;
}
template<>
inline void Z_NR<long>::submul_ui(const Z_NR<long>& a, unsigned long b) {
  data -= a.data * b;
}

template<>
inline void Z_NR<long>::abs(const Z_NR<long>& a) {
  data = labs(a.data);
}

template<>
inline long Z_NR<long>::exponent() const {
  int intExpo;
  double fNorm = frexp(static_cast<double>(data), &intExpo);
  if (data > MAX_LONG_FAST && fabs(fNorm) == 0.5)
    return computeLongExponent(data);
  else
    return static_cast<long>(intExpo);
}

template<>
inline void Z_NR<long>::swap(Z_NR<long>& a) {
  std::swap(data, a.data);
}

#ifdef FPLLL_V3_COMPAT

inline int sizeinbase2(long data) {
  long y = abs(data);
  int resul = 0;
  if (y == 0) resul=1;
  else {
    while (y > 1) {
      resul++;
      y >>= 1;  //y /= 2;
    }
  }
  return resul;
}

template<>
inline void Z_NR<long>::print() const {
  cout << data;
}
template<>
inline void Z_NR<long>::printerr() const {
  cerr << data;
}
template<>
inline void Z_NR<long>::read() {
  cin >> data;
}
template<>
inline double Z_NR<long>::get_d_2exp(long* expo) const {
  int intExpo = 0;
  double x = frexp(static_cast<double>(data), &intExpo);
  *expo = intExpo;
  return x;
}
template<>
inline void Z_NR<long>::set(/*const*/ long& s) {
  data = s;
}
template<>
inline void Z_NR<long>::set(unsigned long s) {
  data = static_cast<long>(s);
}

#endif // #ifdef FPLLL_V3_COMPAT
#endif // #ifdef FPLLL_WITH_ZLONG

/***************************
 * Z=double specialization *
 ***************************/

#ifdef FPLLL_WITH_ZDOUBLE

template<>
inline Z_NR<double>::Z_NR() {}

template<>
inline Z_NR<double>::Z_NR(const Z_NR<double>& z) : data(z.data) {}

template<>
inline Z_NR<double>::~Z_NR() {}

template<>
inline double Z_NR<double>::get_d() const {
  return data;
}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<>
inline long double Z_NR<double>::get_ld() const {
  return static_cast<long double>(data);
}
#endif
template<>
inline long Z_NR<double>::get_si() const {
  return static_cast<long>(rint(data));
}
template<>
inline void Z_NR<double>::set_str(const char* s) {
  data = atof(s);
}
template<>
inline void Z_NR<double>::operator=(const Z_NR<double>& a) {
  data = a.data;
}
template<>
inline void Z_NR<double>::operator=(long a) {
  data = static_cast<double>(a);
}

template<>
inline int Z_NR<double>::sgn() const {
  if (data > 0.0) return 1;
  if (data == 0.0) return 0;
  return -1;
}
template<>
inline int Z_NR<double>::cmp(const Z_NR<double>& m) const {
  if (data > m.data) return 1;
  if (data == m.data) return 0;
  return -1;
}
template<>
inline bool Z_NR<double>::operator<(const Z_NR<double>& a) const {
  return data < a.data;
}
template<>
inline bool Z_NR<double>::operator<(long a) const {
  return data < static_cast<double>(a);
}
template<>
inline bool Z_NR<double>::operator>(const Z_NR<double>& a) const {
  return data > a.data;
}
template<>
inline bool Z_NR<double>::operator>(long a) const {
  return data > static_cast<double>(a);
}
template<>
inline bool Z_NR<double>::operator<=(const Z_NR<double>& a) const {
  return data <= a.data;
}
template<>
inline bool Z_NR<double>::operator<=(long a) const {
  return data <= static_cast<double>(a);
}
template<>
inline bool Z_NR<double>::operator>=(const Z_NR<double>& a) const {
  return data >= a.data;
}
template<>
inline bool Z_NR<double>::operator>=(long a) const {
  return data >= static_cast<double>(a);
}
template<>
inline bool Z_NR<double>::operator==(const Z_NR<double>& a) const {
  return data == a.data;
}
template<>
inline bool Z_NR<double>::operator==(long a) const {
  return data == static_cast<double>(a);
}
template<>
inline bool Z_NR<double>::operator!=(const Z_NR<double>& a) const {
  return data != a.data;
}
template<>
inline bool Z_NR<double>::operator!=(long a) const {
  return data != static_cast<double>(a);
}

template<>
inline void Z_NR<double>::add(const Z_NR<double>& a, const Z_NR<double>& b) {
  data = a.data + b.data;
}
template<>
inline void Z_NR<double>::add_ui(const Z_NR<double>& a, unsigned int b) {
  data = a.data + static_cast<double>(b);
}
template<>
inline void Z_NR<double>::sub(const Z_NR<double>& a, const Z_NR<double>& b) {
  data = a.data - b.data;
}
template<>
inline void Z_NR<double>::neg(const Z_NR<double>& a) {
  data = -a.data;
}
template<>
inline void Z_NR<double>::mul(const Z_NR<double>& a, const Z_NR<double>& b) {
  data = a.data * b.data;
}
template<>
inline void Z_NR<double>::mul_si(const Z_NR<double>& a, long int b) {
  data = a.data * static_cast<double>(b);
}
template<>
inline void Z_NR<double>::mul_ui(const Z_NR<double>& a, long unsigned int b) {
  data = a.data * static_cast<double>(b);
}
template<>
inline void Z_NR<double>::mul_2si(const Z_NR<double>& a, long b) {
  data = ldexp(a.data, b);
}
template<>
inline void Z_NR<double>::div_2si(const Z_NR<double>& a, long b) {
  data = ldexp(a.data, -b);
}
template<>
inline void Z_NR<double>::addmul(const Z_NR<double>& a, const Z_NR<double>& b) {
  data += a.data * b.data;
} 
template<>
inline void Z_NR<double>::addmul_ui(const Z_NR<double>& a, unsigned long b) {
  data += a.data * static_cast<double>(b);
}
template<>
inline void Z_NR<double>::addmul_si(const Z_NR<double>& a, long b) {
  data += a.data * static_cast<double>(b);
}
template<>
inline void Z_NR<double>::submul(const Z_NR<double>& a, const Z_NR<double>& b) {
  data -= a.data * b.data;
} 
template<>
inline void Z_NR<double>::submul_ui(const Z_NR<double>& a, unsigned long b) {
  data -= a.data * static_cast<double>(b);
}

template<>
inline void Z_NR<double>::abs(const Z_NR<double>& a) {
  data = fabs(a.data);
}
template<>
inline long Z_NR<double>::exponent() const {
  int intExpo;
  frexp(data, &intExpo);
  return static_cast<long>(intExpo);
}
template<>
inline void Z_NR<double>::swap(Z_NR<double>& a) {
  std::swap(data, a.data);
}

#ifdef FPLLL_V3_COMPAT

template<>
inline void Z_NR<double>::print() const {
  cout << data;
}
template<>
inline void Z_NR<double>::printerr() const {
  cerr << data;
}
template<>
inline void Z_NR<double>::read() {
  cin >> data;
}
template<>
inline double Z_NR<double>::get_d_2exp(long* expo) const {
  int intExpo = 0;
  double x = frexp(data, &intExpo);
  *expo = intExpo;
  return x;
}
template<>
inline void Z_NR<double>::set(/*const*/ double& s) {
  data = s;
}
template<>
inline void Z_NR<double>::set(unsigned long s) {
  data = static_cast<double>(s);
}

#endif // #ifdef FPLLL_V3_COMPAT
#endif // #ifdef FPLLL_WITH_ZDOUBLE

/**********************
 * mpz specialization *
 **********************/

template<>
inline Z_NR<mpz_t>::Z_NR() {
  mpz_init(data);
}
template<>
inline Z_NR<mpz_t>::Z_NR(const Z_NR<mpz_t>& z) {
  mpz_init_set(data, z.data);
}
template<>
inline Z_NR<mpz_t>::~Z_NR() {
  mpz_clear(data);
}

template<>
inline double Z_NR<mpz_t>::get_d() const {
  return mpz_get_d(data);
}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<>
inline long double Z_NR<mpz_t>::get_ld() const {
  return LDConvHelper::mpz_get_ld(data);
}
#endif
template<>
inline long Z_NR<mpz_t>::get_si() const {
  return mpz_get_si(data);
}
template<>
inline void Z_NR<mpz_t>::set_str(const char* s) {
  mpz_set_str(data, s, 10);
}
template<>
inline void Z_NR<mpz_t>::operator=(const Z_NR<mpz_t>& a) {
  mpz_set(data, a.data);
}
template<>
inline void Z_NR<mpz_t>::operator=(long a) {
  mpz_set_si(data, a);
}

template<>
inline int Z_NR<mpz_t>::sgn() const {
  return mpz_sgn(data);
}
template<>
inline int Z_NR<mpz_t>::cmp(const Z_NR<mpz_t>& m) const {
  int c = mpz_cmp(data, m.data);
  if (c > 0) return 1;
  if (c == 0) return 0;
  return -1;
}
template<>
inline bool Z_NR<mpz_t>::operator<(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) < 0;
}
template<>
inline bool Z_NR<mpz_t>::operator<(long a) const {
  return mpz_cmp_si(data, a) < 0;
}
template<>
inline bool Z_NR<mpz_t>::operator>(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) > 0;
}
template<>
inline bool Z_NR<mpz_t>::operator>(long a) const {
  return mpz_cmp_si(data, a) > 0;
}
template<>
inline bool Z_NR<mpz_t>::operator<=(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) <= 0;
}
template<>
inline bool Z_NR<mpz_t>::operator<=(long a) const {
  return mpz_cmp_si(data, a) <= 0;
}
template<>
inline bool Z_NR<mpz_t>::operator>=(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) >= 0;
}
template<>
inline bool Z_NR<mpz_t>::operator>=(long a) const {
  return mpz_cmp_si(data, a) >= 0;
}
template<>
inline bool Z_NR<mpz_t>::operator==(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) == 0;
}
template<>
inline bool Z_NR<mpz_t>::operator==(long a) const {
  return mpz_cmp_si(data, a) == 0;
}
template<>
inline bool Z_NR<mpz_t>::operator!=(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) != 0;
}
template<>
inline bool Z_NR<mpz_t>::operator!=(long a) const {
  return mpz_cmp_si(data, a) != 0;
}

template<>
inline void Z_NR<mpz_t>::add(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_add(data, a.data, b.data);
}
template<>
inline void Z_NR<mpz_t>::add_ui(const Z_NR<mpz_t>& a, unsigned int b) {
  mpz_add_ui(data, a.data, b);
}
template<>
inline void Z_NR<mpz_t>::sub(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_sub(data, a.data, b.data);
}
template<>
inline void Z_NR<mpz_t>::neg(const Z_NR<mpz_t>& a) {
  mpz_neg(data, a.data);
}
template<>
inline void Z_NR<mpz_t>::mul(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_mul(data, a.data, b.data);
}
template<>
inline void Z_NR<mpz_t>::mul_si(const Z_NR<mpz_t>& a, long b) {
  mpz_mul_si(data, a.data, b);
}
template<>
inline void Z_NR<mpz_t>::mul_ui(const Z_NR<mpz_t>& a, unsigned long b) {
  mpz_mul_ui(data, a.data, b);
}
template<>
inline void Z_NR<mpz_t>::mul_2si(const Z_NR<mpz_t>& a, long b) {
  if (b >= 0)
    mpz_mul_2exp(data, a.data, b);
  else
    mpz_div_2exp(data, a.data, -b);
}
template<>
inline void Z_NR<mpz_t>::div_2si(const Z_NR<mpz_t>& a, long b) {
  if (b >= 0)
    mpz_div_2exp(data, a.data, b);
  else
    mpz_mul_2exp(data, a.data, -b);
}
template<>
inline void Z_NR<mpz_t>::addmul(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_addmul(data, a.data, b.data);
}
template<>
inline void Z_NR<mpz_t>::addmul_ui(const Z_NR<mpz_t>& a, unsigned long b) {
  mpz_addmul_ui(data, a.data, b);
}
template<>
inline void Z_NR<mpz_t>::addmul_si(const Z_NR<mpz_t>& a, long b) {
  if (b >= 0)
    mpz_addmul_ui(data, a.data, static_cast<unsigned long>(b));
  else
    mpz_submul_ui(data, a.data, static_cast<unsigned long>(-b));
}
template<>
inline void Z_NR<mpz_t>::submul(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_submul(data, a.data, b.data);
}
template<>
inline void Z_NR<mpz_t>::submul_ui(const Z_NR<mpz_t>& a,unsigned long b) {
  mpz_submul_ui(data, a.data, b);
}

template<>
inline void Z_NR<mpz_t>::abs(const Z_NR<mpz_t>& x) {
  mpz_abs(data,x.data);
}
template<>
inline long Z_NR<mpz_t>::exponent() const {
  long expo;
  mpz_get_d_2exp(&expo, data);
  return expo;
}
template<>
inline void Z_NR<mpz_t>::randb(int bits) {
  mpz_urandomb(data, RandGen::getGMPState(), bits);
}
template<>
inline void Z_NR<mpz_t>::randm(const Z_NR<mpz_t>& max) {
  mpz_urandomm(data, RandGen::getGMPState(), max.data);
}
template<>
inline void Z_NR<mpz_t>::swap(Z_NR<mpz_t>& a) {
  mpz_swap(data, a.data);
}

#ifdef FPLLL_V3_COMPAT

template<>
inline void Z_NR<mpz_t>::print() const {
  mpz_out_str(stdout, 10, data);
}
template<>
inline void Z_NR<mpz_t>::printerr() const {
  mpz_out_str(stderr, 10, data);
}
template<>
inline void Z_NR<mpz_t>::read() {
  mpz_inp_str(data, stdin, 0);
}
template<>
inline double Z_NR<mpz_t>::get_d_2exp(long* expo) const {
  return mpz_get_d_2exp(expo,data);
}
template<>
inline void Z_NR<mpz_t>::set(/*const*/ mpz_t& d) {
  mpz_set(data,d);
}
template<>
inline void Z_NR<mpz_t>::set(unsigned long d) {
  mpz_set_ui(data,d);
}

#endif // #ifdef FPLLL_V3_COMPAT

/*********
 * FP_NR *
 *********/

template<class F>
inline void FP_NR<F>::addmul(const FP_NR<F>& b, const FP_NR<F>& c, mp_rnd_t rnd) {
  FP_NR<F> product;
  product.mul(b, c, rnd);
  add(*this, product, rnd);
}
template<class F>
inline void FP_NR<F>::submul(const FP_NR<F>& b, const FP_NR<F>& c, mp_rnd_t rnd) {
  FP_NR<F> product;
  product.mul(b, c, rnd);
  sub(*this, product, rnd);
}

/*************************
 * double specialization *
 *************************/

template<>
inline FP_NR<double>::FP_NR() {}
template<>
inline FP_NR<double>::FP_NR(const FP_NR<double>& f) {
  data = f.data;
}
template<>
inline FP_NR<double>::~FP_NR() {}

template<>
inline long FP_NR<double>::exponent() const {
  return ilogb(data) + 1;
}

template<>
inline double FP_NR<double>::get_d(mp_rnd_t /*rnd*/) const {
  return static_cast<double>(data);
}
template<>
inline long FP_NR<double>::get_si() const {
  return static_cast<long>(data);
}
template<>
inline long FP_NR<double>::get_si_exp_we(long& expo, long expoAdd) const {
  if (data == 0)
    expo = 0;
  else
    expo = max(exponent() + expoAdd - numeric_limits<long>::digits, 0L);
  return static_cast<long>(ldexp(data, expoAdd - expo));
}
template<>
inline long FP_NR<double>::get_si_exp(long& expo) const {
  return get_si_exp_we(expo, 0);
}
template<>
inline void FP_NR<double>::operator=(const FP_NR<double>& f) {
  data = f.data;
}
template<>
inline void FP_NR<double>::operator=(double d) {
  data = d;
}

template<>
inline int  FP_NR<double>::cmp(const FP_NR<double>& b) const {
  if (data > b.data) return 1;
  if (data < b.data) return -1;
  return 0;
}
template<>
inline int  FP_NR<double>::cmp(double b) const {
  if (data > b) return 1;
  if (data < b) return -1;
  return 0;
}
template<>
inline int FP_NR<double>::sgn() const {
  if (data > 0) return 1;
  if (data < 0) return -1;
  return 0;
}
template<>
inline bool FP_NR<double>::operator<=(const FP_NR<double>& a) const {
  return data <= a.data;
}
template<>
inline bool FP_NR<double>::operator<=(double a) const {
  return data <= a;
}
template<>
inline bool FP_NR<double>::operator>=(const FP_NR<double>& a) const {
  return data >= a.data;
}
template<>
inline bool FP_NR<double>::operator>=(double a) const {
  return data >= a;
}
template<>
inline bool FP_NR<double>::operator<(const FP_NR<double>& a) const {
  return data < a.data;
}
template<>
inline bool FP_NR<double>::operator<(double a) const {
  return data < a;
}
template<>
inline bool FP_NR<double>::operator>(const FP_NR<double>& a) const {
  return data > a.data;
}
template<>
inline bool FP_NR<double>::operator>(double a) const {
  return data > a;
}

template<>
inline void FP_NR<double>::add(const FP_NR<double>& b, const FP_NR<double>& c, mp_rnd_t /*rnd*/) {
  data = b.data + c.data;
}
template<>
inline void FP_NR<double>::sub(const FP_NR<double>& b, const FP_NR<double>& c, mp_rnd_t /*rnd*/) {
  data = b.data - c.data;
}
template<>
inline void FP_NR<double>::neg(const FP_NR<double>& b) {
  data = -b.data;
}
template<>
inline void FP_NR<double>::mul(const FP_NR<double>& b, const FP_NR<double>& c, mp_rnd_t /*rnd*/) {
  data = b.data * c.data;
}
template<>
inline void FP_NR<double>::mul_2si(const FP_NR<double>& b, long c) {
  data = ldexp(b.data, static_cast<int>(c));
}
template<>
inline void FP_NR<double>::div(const FP_NR<double>& b, const FP_NR<double>& c, mp_rnd_t /*rnd*/) {
  data = b.data / c.data;
}

template<>
inline void FP_NR<double>::abs(const FP_NR<double>& b) {
  data = b.data;
  if (data < 0) data = -data;
}
template<>
inline void FP_NR<double>::rnd(const FP_NR<double>& b) {
  data = rint(b.data);
}
template<>
inline void FP_NR<double>::rnd_we(const FP_NR<double>& b, long expoAdd) {
  // If data == 0.0, exponent() is undefined, but both branches will work
  if (b.exponent() + expoAdd >= numeric_limits<double>::digits)
    data = b.data;
  else
    data = ldexp(::rint(ldexp(b.data, expoAdd)), -expoAdd);
}
template<>
inline void FP_NR<double>::floor(const FP_NR<double>& b) {
  data = ::floor(b.data);
}
template<>
inline bool FP_NR<double>::is_zero() const {
  return data == 0;
}
template<>
inline void FP_NR<double>::set_nan() {
  data = NAN;
}
template<>
inline int FP_NR<double>::is_nan() const {
  return data != data;
}
template<>
inline int FP_NR<double>::is_finite() const {
  return isfinite(data);
}
template<>
inline void FP_NR<double>::sqrt(const FP_NR<double>& s, mp_rnd_t /*rnd*/) {
  data = ::sqrt(s.data);
}
template<>
inline void FP_NR<double>::pow_si(const FP_NR<double>& a, long b, mp_rnd_t /*rnd*/) {
  data = ::pow(a.data, static_cast<double>(b));
}
template<>
inline void FP_NR<double>::exponential(const FP_NR<double>& a, mp_rnd_t /*rnd*/) {
  data = ::exp(a.data);
}
template<>
inline void FP_NR<double>::log(const FP_NR<double>& a, mp_rnd_t /*rnd*/) {
  data = ::log(a.data);
}
template<>
inline void FP_NR<double>::swap(FP_NR<double>& a) {
  std::swap(data, a.data);
}

template<>
inline unsigned int FP_NR<double>::getprec() {
  return PREC_DOUBLE;
}

template<>
inline unsigned int FP_NR<double>::setprec(unsigned int /*prec*/) {
  // ignored
  return getprec();
}

#ifdef FPLLL_V3_COMPAT

template<>
inline void FP_NR<double>::print() const {
  cout << data;
}
template<>
inline void FP_NR<double>::printerr() const {
  cerr << data;
}
template<>
inline void FP_NR<double>::set(unsigned int s) {
  data = static_cast<double>(s);
}

#endif // #ifdef FPLLL_V3_COMPAT


#ifdef FPLLL_WITH_LONG_DOUBLE
/******************************
 * long double specialization *
 ******************************/

template<>
inline FP_NR<long double>::FP_NR() {}
template<>
inline FP_NR<long double>::FP_NR(const FP_NR<long double>& f) : data (f.data) {}
template<>
inline FP_NR<long double>::~FP_NR() {}

template<>
inline long FP_NR<long double>::exponent() const {
  return static_cast<long>(ilogbl(data) + 1);
}

template<>
inline double FP_NR<long double>::get_d(mp_rnd_t /*rnd*/) const {
  return static_cast<double>(data);
}
template<>
inline long FP_NR<long double>::get_si() const {
  return static_cast<long>(data);
}
template<>
inline long FP_NR<long double>::get_si_exp_we(long& expo, long expoAdd) const {
  if (data == 0.0)
    expo = 0;
  else
    expo = max(exponent() + expoAdd - numeric_limits<long>::digits, 0L);
  return static_cast<long>(ldexpl(data, expoAdd - expo));
}
template<>
inline long FP_NR<long double>::get_si_exp(long& expo) const {
  return get_si_exp_we(expo, 0);
}
template<>
inline void FP_NR<long double>::operator=(const FP_NR<long double>& f) {
  data = f.data;
}
template<>
inline void FP_NR<long double>::operator=(double d) {
  data = d;
}

template<>
inline int FP_NR<long double>::cmp(const FP_NR<long double>& b) const {
  if (data > b.data) return 1;
  if (data < b.data) return -1;
  return 0;
}
template<>
inline int FP_NR<long double>::cmp(double b) const {
  if (data > b) return 1;
  if (data < b) return -1;
  return 0;
}
template<>
inline int FP_NR<long double>::sgn() const {
  if (data > 0) return 1;
  if (data < 0) return -1;
  return 0;
}
template<>
inline bool FP_NR<long double>::operator<=(const FP_NR<long double>& a) const {
  return data <= a.data;
}
template<>
inline bool FP_NR<long double>::operator<=(double a) const {
  return data <= a;
}
template<>
inline bool FP_NR<long double>::operator>=(const FP_NR<long double>& a) const {
  return data >= a.data;
}
template<>
inline bool FP_NR<long double>::operator>=(double a) const {
  return data >= a;
}
template<>
inline bool FP_NR<long double>::operator<(const FP_NR<long double>& a) const {
  return data < a.data;
}
template<>
inline bool FP_NR<long double>::operator<(double a) const {
  return data < a;
}
template<>
inline bool FP_NR<long double>::operator>(const FP_NR<long double>& a) const {
  return data > a.data;
}
template<>
inline bool FP_NR<long double>::operator>(double a) const {
  return data > a;
}

template<>
inline void FP_NR<long double>::add(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t /*rnd*/) {
  data = b.data + c.data;
}
template<>
inline void FP_NR<long double>::sub(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t /*rnd*/) {
  data = b.data - c.data;
}
template<>
inline void FP_NR<long double>::neg(const FP_NR<long double>& b) {
  data = -b.data;
}
template<>
inline void FP_NR<long double>::mul(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t /*rnd*/) {
  data = b.data * c.data;
}
template<>
inline void FP_NR<long double>::mul_2si(const FP_NR<long double>& b, long c) {
  data = ldexpl(b.data, static_cast<int>(c));
}
template<>
inline void FP_NR<long double>::div(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t /*rnd*/) {
  data = b.data / c.data;
}

template<>
inline void FP_NR<long double>::abs(const FP_NR<long double>& b) {
  data = b.data;
  if (data < 0) data = -data;
}
template<>
inline void FP_NR<long double>::rnd(const FP_NR<long double>& b) {
  data = rintl(b.data);
}
template<>
inline void FP_NR<long double>::rnd_we(const FP_NR<long double>& b, long expoAdd) {
  // If data == 0.0, exponent() is undefined, but both branches will work
  if (b.exponent() + expoAdd >= numeric_limits<long double>::digits)
    data = b.data;
  else
    data = ldexpl(rintl(ldexpl(b.data, expoAdd)), -expoAdd);
}
template<>
inline void FP_NR<long double>::floor(const FP_NR<long double>& b) {
  data = floorl(b.data);
}

template<>
inline bool FP_NR<long double>::is_zero() const {
  return data == 0;
}
template<>
inline void FP_NR<long double>::set_nan() {
  data = NAN;
}
template<>
inline int FP_NR<long double>::is_nan() const {
  return data != data;
}
template<>
inline int FP_NR<long double>::is_finite() const {
  return isfinite(data);
}
template<>
inline void FP_NR<long double>::sqrt(const FP_NR<long double>& s, mp_rnd_t /*rnd*/) {
  data = sqrtl(s.data);
}
template<>
inline void FP_NR<long double>::pow_si(const FP_NR<long double>& a, long b, mp_rnd_t /*rnd*/) {
  data = powl(a.data, static_cast<long double>(b));
}
template<>
inline void FP_NR<long double>::exponential(const FP_NR<long double>& a, mp_rnd_t /*rnd*/) {
  data = expl(a.data);
}
template<>
inline void FP_NR<long double>::log(const FP_NR<long double>& a, mp_rnd_t /*rnd*/) {
  data = logl(a.data);
}
template<>
inline void FP_NR<long double>::swap(FP_NR<long double>& a) {
  std::swap(data, a.data);
}

template<>
inline unsigned int FP_NR<long double>::getprec() {
  return numeric_limits<long double>::digits;
}
template<>
inline unsigned int FP_NR<long double>::setprec(unsigned int /*prec*/) {
  return getprec(); // ignored
}

#ifdef FPLLL_V3_COMPAT

template<>
inline void FP_NR<long double>::print() const {
  cout << data;
}
template<>
inline void FP_NR<long double>::printerr() const {
  cerr << data;
}
template<>
inline void FP_NR<long double>::set(unsigned int s) {
  data = static_cast<long double>(s);
}
#endif // #ifdef FPLLL_V3_COMPAT

#endif // #ifef FPLLL_WITH_LONG_DOUBLE


/*************
 * DPE spec. *
 *************/

#ifdef FPLLL_WITH_DPE

template<>
inline FP_NR<dpe_t>::FP_NR() {
  dpe_init(data);
}
template<>
inline FP_NR<dpe_t>::FP_NR(const FP_NR<dpe_t>& f) {
  dpe_init(data);
  dpe_set(data, f.data);
}
template<>
inline FP_NR<dpe_t>::~FP_NR() {
  dpe_clear(data);
}

template<>
inline double FP_NR<dpe_t>::get_d(mp_rnd_t /*rnd*/) const {
  return dpe_get_d(data);
}
template<>
inline long FP_NR<dpe_t>::get_si() const {
  return dpe_get_si(data);
}
template<>
inline long FP_NR<dpe_t>::get_si_exp(long& expo) const {
  long result;
  expo = dpe_get_si_exp(&result, data);
  if (dpe_zero_p(data))
    expo = 0;
  else if (expo < 0) {
    /* NOTE: conversion of result to double is exact even if
        sizeof(long) = 8 */
    result = static_cast<long>(ldexp(static_cast<double>(result), expo));
    expo = 0;
  }
  return result;
}
template<>
inline long FP_NR<dpe_t>::get_si_exp_we(long& expo, long /*expoAdd*/) const {
  return get_si_exp(expo); // NOTE: expoAdd = 0
}
template<>
inline void FP_NR<dpe_t>::operator=(const FP_NR<dpe_t>& f) {
  dpe_set(data, f.data);
}
template<>
inline void FP_NR<dpe_t>::operator=(double d) {
  dpe_set_d(data, d);
}

template<>
inline void FP_NR<dpe_t>::log(const FP_NR<dpe_t>& a, mp_rnd_t /*rnd*/) {
  dpe_ugly_log(data, a.data);
  /*  FPLLL_DEBUG_ABORT("FP_NR<dpe_t>::log not implemented"); */
}

template<>
inline void FP_NR<dpe_t>::exponential(const FP_NR<dpe_t>& a, mp_rnd_t /*rnd*/) {
  dpe_ugly_exp(data, a.data);
  /*  FPLLL_DEBUG_ABORT("FP_NR<dpe_t>::exp not implemented"); */
}


template<>
inline int FP_NR<dpe_t>::cmp(const FP_NR<dpe_t>& a) const {
  return dpe_cmp(data, a.data);
}
template<>
inline int FP_NR<dpe_t>::cmp(double a) const {
  return dpe_cmp_d(data, a);
}
template<>
inline int FP_NR<dpe_t>::sgn() const {
  return cmp(0.0);
}
template<>
inline bool FP_NR<dpe_t>::operator<=(const FP_NR<dpe_t>& a) const {
  return dpe_cmp(data, a.data) <= 0;
}
template<>
inline bool FP_NR<dpe_t>::operator<=(double a) const {
  return dpe_cmp_d(data, a) <= 0;
}
template<>
inline bool FP_NR<dpe_t>::operator>=(const FP_NR<dpe_t>& a) const {
  return dpe_cmp(data, a.data) >= 0;
}
template<>
inline bool FP_NR<dpe_t>::operator>=(double a) const {
  return dpe_cmp_d(data, a) >= 0;
}

template<>
inline void FP_NR<dpe_t>::add(const FP_NR<dpe_t>& a, const FP_NR<dpe_t>& b, mp_rnd_t /*rnd*/) {
  dpe_add(data, a.data, b.data);
}
template<>
inline void FP_NR<dpe_t>::sub(const FP_NR<dpe_t>& a, const FP_NR<dpe_t>& b, mp_rnd_t /*rnd*/) {
  dpe_sub(data, a.data, b.data);
}
template<>
inline void FP_NR<dpe_t>::neg(const FP_NR<dpe_t>& a) {
  dpe_neg(data, a.data);
}
template<>
inline void FP_NR<dpe_t>::mul(const FP_NR<dpe_t>& a, const FP_NR<dpe_t>& b, mp_rnd_t /*rnd*/) {
  dpe_mul(data, a.data, b.data);
}
template<>
inline void FP_NR<dpe_t>::mul_2si(const FP_NR<dpe_t>& a, long b) {
  dpe_mul_2si(data, a.data, b);
}
template<>
inline void FP_NR<dpe_t>::div(const FP_NR<dpe_t>& a, const FP_NR<dpe_t>& b, mp_rnd_t /*rnd*/) {
  dpe_div(data, a.data, b.data);
}

template<>
inline void FP_NR<dpe_t>::abs(const FP_NR<dpe_t>& a) {
  dpe_abs(data, a.data);
}
template<>
inline void FP_NR<dpe_t>::rnd(const FP_NR<dpe_t>& a) {
  dpe_round(data, a.data);
}
template<>
inline void FP_NR<dpe_t>::rnd_we(const FP_NR<dpe_t>& a, long /*expoAdd*/) {
  dpe_round(data, a.data); // NOTE: expoAdd = 0
}
template<>
inline void FP_NR<dpe_t>::floor(const FP_NR<dpe_t>& a) {
  dpe_floor(data, a.data);
}
template<>
inline long FP_NR<dpe_t>::exponent() const {
  return DPE_EXP(data);
}
template<>
inline bool FP_NR<dpe_t>::is_zero() const {
  return dpe_zero_p(data);
}
template<>
inline void FP_NR<dpe_t>::set_nan() {
  //dpe_set_d(data, NAN); // DPE_UNLIKELY branch in dpe_normalize
  DPE_MANT(data) = NAN;
}
template<>
inline int FP_NR<dpe_t>::is_nan() const {
  return DPE_MANT(data) != DPE_MANT(data);
}
template<>
inline int FP_NR<dpe_t>::is_finite() const {
  return isfinite(DPE_MANT(data));
}
template<>
inline void FP_NR<dpe_t>::sqrt(const FP_NR<dpe_t>& a, mp_rnd_t /*rnd*/) {
  dpe_sqrt(data, a.data);
}
template<>
inline void FP_NR<dpe_t>::swap(FP_NR<dpe_t>& a) {
  dpe_swap(data, a.data);
}

template<>
inline unsigned int FP_NR<dpe_t>::getprec() {
  return DPE_BITSIZE;
}
template<>
inline unsigned int FP_NR<dpe_t>::setprec(unsigned int /*prec*/) {
  return getprec(); // ignored
}

#ifdef FPLLL_V3_COMPAT

template<>
inline void FP_NR<dpe_t>::print() const {
  dpe_out_str(stdout,10, data);
  fflush(stdout);
}
template<>
inline void FP_NR<dpe_t>::printerr() const {
  dpe_out_str(stderr,10, data);
  fflush(stderr);
}
template<>
inline void FP_NR<dpe_t>::set(unsigned int s) {
  dpe_set_d(data, static_cast<double>(s));
}

#endif // #ifdef FPLLL_V3_COMPAT
#endif // #ifdef FPLLL_WITH_DPE

/*************
 * MPFR spec. *
 *************/

template<>
inline FP_NR<mpfr_t>::FP_NR() {
  mpfr_init(data);
}
template<>
inline FP_NR<mpfr_t>::FP_NR(const FP_NR<mpfr_t>& f) {
  mpfr_init_set(data, f.data, GMP_RNDN);
}
template<>
inline FP_NR<mpfr_t>::~FP_NR() {
  mpfr_clear(data);
}

template<>
inline long FP_NR<mpfr_t>::exponent() const {
  return mpfr_get_exp(data);
}

template<>
inline double FP_NR<mpfr_t>::get_d(mp_rnd_t rnd) const {
  return mpfr_get_d(data, rnd);
}
template<>
inline long FP_NR<mpfr_t>::get_si() const {
  return mpfr_get_si(data,GMP_RNDN);
}
template<>
inline long FP_NR<mpfr_t>::get_si_exp(long& expo) const {
  if (mpfr_zero_p(data)) {
    expo = 0;
  }
  else {
    expo = max(exponent() - numeric_limits<long>::digits, 0L);
  }
  mpfr_t& ncData = const_cast<mpfr_t&>(data);
  mpfr_div_2si(ncData, ncData, expo, GMP_RNDN);
  long result = mpfr_get_si(ncData, GMP_RNDZ);
  mpfr_mul_2si(ncData, ncData, expo, GMP_RNDN);
  return result;
}
template<>
inline long FP_NR<mpfr_t>::get_si_exp_we(long& expo, long /*expoAdd*/) const {
  return get_si_exp(expo);  // NOTE: expoAdd = 0
}

template<>
inline void FP_NR<mpfr_t>::operator=(const FP_NR<mpfr_t>& a) {
  mpfr_set(data, a.data, GMP_RNDN);
}
template<>
inline void FP_NR<mpfr_t>::operator=(double a) {
  mpfr_set_d(data, a, GMP_RNDN);
}

template<>
inline int FP_NR<mpfr_t>::cmp(const FP_NR<mpfr_t>& a) const {
  return mpfr_cmp(data, a.data);
}
template<>
inline int FP_NR<mpfr_t>::cmp(double a) const {
  return mpfr_cmp_d(data, a);
}
template<>
inline int FP_NR<mpfr_t>::sgn() const {
  return mpfr_sgn(data);
}
template<>
inline bool FP_NR<mpfr_t>::operator<=(const FP_NR<mpfr_t>& a) const {
  return mpfr_cmp(data, a.data) <= 0;
}
template<>
inline bool FP_NR<mpfr_t>::operator<=(double a) const {
  return mpfr_cmp_d(data, a) <= 0;
}
template<>
inline bool FP_NR<mpfr_t>::operator>=(const FP_NR<mpfr_t>& a) const {
  return mpfr_cmp(data, a.data) >= 0;
}
template<>
inline bool FP_NR<mpfr_t>::operator>=(double a) const {
  return mpfr_cmp_d(data, a) >= 0;
}

template<>
inline void FP_NR<mpfr_t>::add(const FP_NR<mpfr_t>& a, const FP_NR<mpfr_t>& b, mp_rnd_t rnd) {
  mpfr_add(data, a.data, b.data, rnd);
}
template<>
inline void FP_NR<mpfr_t>::sub(const FP_NR<mpfr_t>& a, const FP_NR<mpfr_t>& b, mp_rnd_t rnd) {
  mpfr_sub(data, a.data, b.data, rnd);
}
template<>
inline void FP_NR<mpfr_t>::neg(const FP_NR<mpfr_t>& a) {
  mpfr_neg(data, a.data, GMP_RNDN);
}
template<>
inline void FP_NR<mpfr_t>::mul(const FP_NR<mpfr_t>& a, const FP_NR<mpfr_t>& b, mp_rnd_t rnd) {
  mpfr_mul(data, a.data, b.data, rnd);
}
template<>
inline void FP_NR<mpfr_t>::mul_2si(const FP_NR<mpfr_t>& a, long b) {
  mpfr_mul_2si(data, a.data, b, GMP_RNDN);
}
template<>
inline void FP_NR<mpfr_t>::div(const FP_NR<mpfr_t>& a, const FP_NR<mpfr_t>& b, mp_rnd_t rnd) {
  mpfr_div(data, a.data, b.data, rnd);
}
template<>
inline void FP_NR<mpfr_t>::addmul(const FP_NR<mpfr_t>& b, const FP_NR<mpfr_t>& c, mp_rnd_t rnd) {
  mpfr_fma(data, b.data, c.data, data, rnd);
}
template<>
inline void FP_NR<mpfr_t>::submul(const FP_NR<mpfr_t>& b, const FP_NR<mpfr_t>& c, mp_rnd_t rnd) {
  mpfr_fms(data, b.data, c.data, data, rnd);
  mpfr_neg(data, data, GMP_RNDN); // Exact
}

template<>
inline void FP_NR<mpfr_t>::abs(const FP_NR<mpfr_t>& a) {
  mpfr_abs(data, a.data, GMP_RNDN);
}
template<>
inline void FP_NR<mpfr_t>::rnd(const FP_NR<mpfr_t>& a) {
  mpfr_round(data, a.data);
}
template<>
inline void FP_NR<mpfr_t>::rnd_we(const FP_NR<mpfr_t>& a, long /*expoAdd*/) {
  mpfr_round(data, a.data); // NOTE: expoAdd = 0
}
template<>
inline void FP_NR<mpfr_t>::floor(const FP_NR<mpfr_t>& a) {
  mpfr_floor(data, a.data);
}
template<>
inline bool FP_NR<mpfr_t>::is_zero() const {
  return mpfr_zero_p(data);
}
template<>
inline void FP_NR<mpfr_t>::set_nan() {
  mpfr_set_nan(data);
}
template<>
inline int FP_NR<mpfr_t>::is_nan() const {
  return mpfr_nan_p(data);
}
template<>
inline int FP_NR<mpfr_t>::is_finite() const {
  return mpfr_number_p(data);
}
template<>
inline void FP_NR<mpfr_t>::sqrt(const FP_NR<mpfr_t>& a, mp_rnd_t rnd) {
  mpfr_sqrt(data, a.data, rnd);
}
template<>
inline void FP_NR<mpfr_t>::pow_si(const FP_NR<mpfr_t>& a, long b, mp_rnd_t rnd) {
  mpfr_pow_si(data, a.data, b, rnd);
}
template<>
inline void FP_NR<mpfr_t>::exponential(const FP_NR<mpfr_t>& a, mp_rnd_t rnd) {
  mpfr_exp(data, a.data, rnd);
}
template<>
inline void FP_NR<mpfr_t>::log(const FP_NR<mpfr_t>& a, mp_rnd_t rnd) {
  mpfr_log(data, a.data, rnd);
}
template<>
inline void FP_NR<mpfr_t>::swap(FP_NR<mpfr_t>& a) {
  mpfr_swap(data, a.data);
}

template<>
inline unsigned int FP_NR<mpfr_t>::getprec() {
  return mpfr_get_default_prec();
}
template<>
inline unsigned int FP_NR<mpfr_t>::setprec(unsigned int prec) {
  int oldprec = getprec();
  mpfr_set_default_prec(prec);
  return oldprec;
}

#ifdef FPLLL_V3_COMPAT

template<>
inline void FP_NR<mpfr_t>::print() const {
  mpfr_out_str(stdout,10,10,data,GMP_RNDN);
  fflush(stdout);
}
template<>
inline void FP_NR<mpfr_t>::printerr() const {
  mpfr_out_str(stderr,10,5,data,GMP_RNDN);
  fflush(stderr);
}
template<>
inline void FP_NR<mpfr_t>::set(unsigned int s) {
  mpfr_set_ui(data,s,GMP_RNDN);
}

#endif // #ifdef FPLLL_V3_COMPAT

template<class F>
inline bool FP_NR<F>::operator<(const FP_NR<F>& a) const {
  return cmp(a) < 0;
}
template<class F>
inline bool FP_NR<F>::operator<(double a) const {
  return cmp(a) < 0;
}
template<class F>
inline bool FP_NR<F>::operator>(const FP_NR<F>& a) const {
  return cmp(a) > 0;
}
template<class F>
inline bool FP_NR<F>::operator>(double a) const {
  return cmp(a) > 0;
}

/* Conversions
   =========== */

/* get_f_exp */

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<double>& f, long& expo) {
  int intExpo;
  f.getData() = frexp(static_cast<double>(data), &intExpo);
  expo = intExpo;
}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<long double>& f, long& expo) {
  int intExpo;
  f.getData() = frexpl(static_cast<long double>(data), &intExpo);
  expo = intExpo;
}
#endif
#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<dpe_t>& /*f*/, long& /*expo*/) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for dpe_t");
}
#endif
template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<mpfr_t>& /*f*/, long& /*expo*/) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for mpfr_t");
}
#endif // #ifdef FPLLL_WITH_ZLONG

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<double>& f, long& expo) {
  int intExpo;
  f.getData() = frexp(data, &intExpo);
  expo = intExpo;
}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<long double>& f, long& expo) {
  int intExpo;
  f.getData() = static_cast<long double>(frexp(data, &intExpo));
  expo = intExpo;
}
#endif
#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<dpe_t>& /*f*/, long& /*expo*/) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for dpe_t");
}
#endif
template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<mpfr_t>& /*f*/, long& /*expo*/) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for mpfr_t");
}
#endif // #ifdef FPLLL_WITH_ZDOUBLE

template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<double>& f, long& expo) {
  f.getData() = mpz_get_d_2exp(&expo, data);
}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<long double>& f, long& expo) {
  f.getData() = LDConvHelper::mpz_get_ld_2exp(&expo, data);
}
#endif
#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<dpe_t>& /*f*/, long& /*expo*/) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for dpe_t");
}
#endif

template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<mpfr_t>& /*f*/, long& /*expo*/) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for mpfr_t");
}

/* set_f */

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void Z_NR<long>::set_f(const FP_NR<double>& a) {
  data = a.get_si();
}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<long>::set_f(const FP_NR<long double>& a) {
  data = a.get_si();
}
#endif
#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<long>::set_f(const FP_NR<dpe_t>& a) {
  data = a.get_si();
}
#endif
template<> template<>
inline void Z_NR<long>::set_f(const FP_NR<mpfr_t>& a) {
  data = a.get_si();
}
#endif // #ifdef FPLLL_WITH_ZLONG

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void Z_NR<double>::set_f(const FP_NR<double>& a) {
  data = a.get_d();
}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<double>::set_f(const FP_NR<long double>& a) {
  data = a.get_d();
}
#endif
#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<double>::set_f(const FP_NR<dpe_t>& a) {
  data = a.get_d();
}
#endif
template<> template<>
inline void Z_NR<double>::set_f(const FP_NR<mpfr_t>& a) {
  data = a.get_d();
}
#endif // #ifdef FPLLL_WITH_ZDOUBLE

template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<double>& a) {
  mpz_set_d(data, a.getData());
}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<long double>& a) {
  LDConvHelper::mpz_set_ld(data, a.getData());
}
#endif
#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<dpe_t>& a) {
  dpe_get_z(data, const_cast<dpe_t&>(a.getData()));
}
#endif
template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<mpfr_t>& a) {
  mpfr_get_z(data, a.getData(), GMP_RNDN);
}

/* get_z_exp */

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<double>::get_z_exp_we(Z_NR<long>& a, long& expo, long expoAdd) const {
  expo = 0;
  a = static_cast<long>(ldexp(data, expoAdd));
}
#endif
#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<double>::get_z_exp_we(Z_NR<double>& a, long& expo, long expoAdd) const {
  expo = 0;
  a.getData() = trunc(ldexp(data, expoAdd));
}
#endif
template<> template<>
inline void FP_NR<double>::get_z_exp_we(Z_NR<mpz_t>& a, long& expo, long expoAdd) const {
  expo = max(exponent() + expoAdd - numeric_limits<double>::digits, 0L);
  mpz_set_d(a.getData(), ldexp(data, expoAdd - expo));
}
template<> template<class Z>
inline void FP_NR<double>::get_z_exp(Z_NR<Z>& a, long& expo) const {
  return get_z_exp_we(a, expo, 0);
}

#ifdef FPLLL_WITH_LONG_DOUBLE

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<long double>::get_z_exp_we(Z_NR<long>& a, long& expo, long expoAdd) const {
  expo = 0;
  a = static_cast<long>(ldexpl(data, expoAdd));
}
#endif
#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<long double>::get_z_exp_we(Z_NR<double>& a, long& expo, long expoAdd) const {
  expo = 0;
  a.getData() = trunc(static_cast<double>(ldexpl(data, expoAdd)));
}
#endif
template<> template<>
inline void FP_NR<long double>::get_z_exp_we(Z_NR<mpz_t>& a, long& expo, long expoAdd) const {
  expo = max(exponent() + expoAdd - numeric_limits<long double>::digits, 0L);
  /* If expo > 0, then
     expoAdd - expo = numeric_limits<long double>::digits - exponent()
     which implies that ldexpl(data, expoAdd - expo) is an integer */
  LDConvHelper::mpz_set_ld(a.getData(), trunc(ldexpl(data, expoAdd - expo)));
}
template<> template<class Z>
inline void FP_NR<long double>::get_z_exp(Z_NR<Z>& a, long& expo) const {
  return get_z_exp_we(a, expo, 0);
}

#endif // #ifdef FPLLL_WITH_LONG_DOUBLE

#ifdef FPLLL_WITH_DPE
#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<dpe_t>::get_z_exp(Z_NR<long>& a, long& expo) const {
  expo = 0;
  a = get_si();
}
#endif
#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<dpe_t>::get_z_exp(Z_NR<double>& a, long& expo) const {
  expo = max(DPE_EXP(data) - numeric_limits<double>::digits, 0);
  a.getData() = trunc(ldexp(DPE_MANT(data), DPE_EXP(data) - expo));
}
#endif
template<> template<>
inline void FP_NR<dpe_t>::get_z_exp(Z_NR<mpz_t>& a, long& expo) const {
  expo = max(DPE_EXP(data) - DPE_BITSIZE, 0);
  mpz_set_d(a.getData(), trunc(ldexp(DPE_MANT(data), DPE_EXP(data) - expo)));
}
template<> template<class Z>
inline void FP_NR<dpe_t>::get_z_exp_we(Z_NR<Z>& a, long& expo, long expoAdd) const {
  return get_z_exp(a, expo);
}
#endif // #ifdef FPLLL_WITH_DPE

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<mpfr_t>::get_z_exp(Z_NR<long>& a, long& expo) const {
  expo = 0;
  a = get_si();
}
#endif
#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<mpfr_t>::get_z_exp(Z_NR<double>& a, long& expo) const {
  expo = 0;
  a.getData() = trunc(mpfr_get_d(data, GMP_RNDZ));
}
#endif
template<> template<>
inline void FP_NR<mpfr_t>::get_z_exp(Z_NR<mpz_t>& a, long& expo) const {
  expo = mpfr_get_z_exp(a.getData(), data);
  if (expo < 0) {
    a.mul_2si(a, expo);
    expo = 0;
  }
}
template<> template<class Z>
inline void FP_NR<mpfr_t>::get_z_exp_we(Z_NR<Z>& a, long& expo, long expoAdd) const {
  return get_z_exp(a, expo);
}

/* set_z */

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<double>::set_z(const Z_NR<long>& a, mp_rnd_t /*rnd*/) {
  data = a.get_d();
}
#endif
#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<double>::set_z(const Z_NR<double>& a, mp_rnd_t /*rnd*/) {
  data = a.get_d();
}
#endif
template<> template<>
inline void FP_NR<double>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t /*rnd*/) {
  data = a.get_d();
}

#ifdef FPLLL_WITH_LONG_DOUBLE
#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<long double>::set_z(const Z_NR<long>& a, mp_rnd_t /*rnd*/) {
  data = a.get_ld();
}
#endif
#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<long double>::set_z(const Z_NR<double>& a, mp_rnd_t /*rnd*/) {
  data = a.get_ld();
}
#endif
template<> template<>
inline void FP_NR<long double>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t /*rnd*/) {
  data = a.get_ld();
}
#endif // #ifdef FPLLL_WITH_LONG_DOUBLE

#ifdef FPLLL_WITH_DPE
#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<dpe_t>::set_z(const Z_NR<long>& a, mp_rnd_t /*rnd*/) {
  dpe_set_si(data, a.getData());
}
#endif
#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<dpe_t>::set_z(const Z_NR<double>& a, mp_rnd_t /*rnd*/) {
  dpe_set_d(data, a.getData());
}
#endif
template<> template<>
inline void FP_NR<dpe_t>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t /*rnd*/) {
  dpe_set_z(data, const_cast<mpz_t&>(a.getData()));
}
#endif // #ifdef FPLLL_WITH_DPE

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<mpfr_t>::set_z(const Z_NR<long>& a, mp_rnd_t rnd) {
  mpfr_set_si(data, a.getData(), rnd);
}
#endif
#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<mpfr_t>::set_z(const Z_NR<double>& a, mp_rnd_t rnd) {
  mpfr_set_d(data, a.getData(), rnd);
}
#endif
template<> template<>
inline void FP_NR<mpfr_t>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t rnd) {
  mpfr_set_z(data, a.getData(), rnd);
}

FPLLL_END_NAMESPACE

#endif
