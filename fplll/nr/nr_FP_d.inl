/****************************
 *  F=double specialization
 ****************************/

#ifndef FPLLL_NR_FP_D_H
#define FPLLL_NR_FP_D_H

#include "../defs.h"
#include "nr_FP.inl"

FPLLL_BEGIN_NAMESPACE


/* Double specialization */


/* constructor */
template<>
inline FP_NR<double>::FP_NR() {}

template<> 
inline FP_NR<double>::FP_NR(const FP_NR<double>& f) {data = f.data;}

template<>
inline FP_NR<double>::~FP_NR() {}

template<>
inline unsigned int FP_NR<double>::get_prec() {return PREC_DOUBLE;}

template<>
inline unsigned int FP_NR<double>::set_prec(unsigned int) {
  return get_prec(); // ignored
}

/* return data */
template<>
inline double FP_NR<double>::get_d(mp_rnd_t) const {
  return static_cast<double>(data);
}

template<>
inline void FP_NR<double>::get_mpfr(mpfr_t r, mp_rnd_t rnd) const {
  mpfr_set_d (r, data, rnd);
}

template<>
inline void FP_NR<double>::set_mpfr(mpfr_t r, mp_rnd_t rnd) {
  data = mpfr_get_d (r, rnd);
}

template<>
inline long FP_NR<double>::get_si() const {
  return static_cast<long>(data);
}

template<>
inline long FP_NR<double>::exponent() const {
  return ilogb(data) + 1;
}

template<>
inline long FP_NR<double>::get_si_exp_we(long& expo, long expo_add) const {
  if (data == 0)
    expo = 0;
  else
    expo = max(exponent() + expo_add - numeric_limits<long>::digits, 0L);
  return static_cast<long>(ldexp(data, expo_add - expo));
}

template<>
inline long FP_NR<double>::get_si_exp(long& expo) const {
  return get_si_exp_we(expo, 0);
}

/*  comparison */
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

/* operators */
template<>
inline void FP_NR<double>::operator=(const FP_NR<double>& f) {
  data = f.data;
}

template<>
inline void FP_NR<double>::operator=(double d) {
  data = d;
}

template<>
inline void FP_NR<double>::operator=(const char *s) {
  data = strtod(s, NULL);
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
inline bool FP_NR<double>::is_zero() const {
  return data == 0;
}

template<>
inline int FP_NR<double>::is_nan() const {
  return data != data;
}

template<>
inline int FP_NR<double>::is_finite() const {
  return isfinite(data);
}

/* arithmetic */
template<>
inline void FP_NR<double>::add(const FP_NR<double>& b, const FP_NR<double>& c, mp_rnd_t /*rnd*/) {
  data = b.data + c.data;
}

template<>
inline void FP_NR<double>::sub(const FP_NR<double>& b, const FP_NR<double>& c, mp_rnd_t /*rnd*/) {
  data = b.data - c.data;
}

template<>
inline void FP_NR<double>::mul(const FP_NR<double>& b, const FP_NR<double>& c, mp_rnd_t /*rnd*/) {
  data = b.data * c.data;
}

template<>
inline void FP_NR<double>::mul_d(const FP_NR<double>& b, const double c, mp_rnd_t /*rnd*/) {
  data = b.data * c;
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
inline void FP_NR<double>::addmul(const FP_NR<double>& b, const FP_NR<double>& c, mp_rnd_t /*rnd*/) {
  data = data + b.data * c.data;
}

template<>
inline void FP_NR<double>::submul(const FP_NR<double>& b, const FP_NR<double>& c, mp_rnd_t /*rnd*/) {
  data = data - b.data * c.data;
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
inline void FP_NR<double>::sqrt(const FP_NR<double>& s, mp_rnd_t /*rnd*/) {
  data = ::sqrt(s.data);
}

template<>
inline void FP_NR<double>::root(const FP_NR<double>& a, unsigned int k, mp_rnd_t /*rnd*/) {
  data = ::pow(a.data, 1/(static_cast<double>(k)));
}

template<>
inline void FP_NR<double>::neg(const FP_NR<double>& b) {
  data = -b.data;
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
inline void FP_NR<double>::rnd_we(const FP_NR<double>& b, long expo_add) {
  // If data == 0.0, exponent() is undefined, but both branches will work
  if (b.exponent() + expo_add >= numeric_limits<double>::digits)
    data = b.data;
  else
    data = ldexp(::rint(ldexp(b.data, expo_add)), -expo_add);
}

template<>
inline void FP_NR<double>::floor(const FP_NR<double>& b) {
  data = ::floor(b.data);
}

template<>
inline void FP_NR<double>::set_nan() {
  data = NAN;
}

template<>
inline void FP_NR<double>::swap(FP_NR<double>& a) {
  std::swap(data, a.data);
}


FPLLL_END_NAMESPACE

#endif
