/****************************
 *  F=long double specialization
 ****************************/

#ifndef FPLLL_NR_FP_LD_H
#define FPLLL_NR_FP_LD_H

#include "../defs.h"
#include "nr_FP.inl"

FPLLL_BEGIN_NAMESPACE

/* LD specialization if defined */
#ifdef FPLLL_WITH_LONG_DOUBLE


/* constructor */
template<>
inline FP_NR<long double>::FP_NR() {}

template<>
inline FP_NR<long double>::FP_NR(const FP_NR<long double>& f) : data (f.data) {}

template<>
inline FP_NR<long double>::~FP_NR() {}

template<>
inline unsigned int FP_NR<long double>::get_prec() {
  return numeric_limits<long double>::digits;
}

template<>
inline unsigned int FP_NR<long double>::set_prec(unsigned int /*prec*/) {
  return get_prec(); // ignored
}

/* return data */
template<>
inline double FP_NR<long double>::get_d(mp_rnd_t /*rnd*/) const {
  return static_cast<double>(data);
}

template<>
inline void FP_NR<long double>::get_mpfr(mpfr_t r, mp_rnd_t rnd) const {
  mpfr_set_ld (r, data, rnd);
}

template<>
inline void FP_NR<long double>::set_mpfr(mpfr_t r, mp_rnd_t rnd) {
  data = mpfr_get_ld (r, rnd);
}

template<>
inline long FP_NR<long double>::get_si() const {
  return static_cast<long>(data);
}

template<>
inline long FP_NR<long double>::exponent() const {
  return static_cast<long>(ilogbl(data) + 1);
}

template<>
inline long FP_NR<long double>::get_si_exp_we(long& expo, long expo_add) const {
  if (data == 0.0)
    expo = 0;
  else
    expo = max(exponent() + expo_add - numeric_limits<long>::digits, 0L);
  return static_cast<long>(ldexpl(data, expo_add - expo));
}

template<>
inline long FP_NR<long double>::get_si_exp(long& expo) const {
  return get_si_exp_we(expo, 0);
}

/*  comparison */
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

/* operators */
template<>
inline FP_NR<long double>& FP_NR<long double>::operator=(const FP_NR<long double>& f) {
  data = f.data;
  return *this;
}

template<>
inline FP_NR<long double>& FP_NR<long double>::operator=(double d) {
  data = d;
  return *this;
}

template<>
inline FP_NR<long double>& FP_NR<long double>::operator=(const char *s) {
  data = strtold(s, NULL);
  return *this;
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
inline bool FP_NR<long double>::is_zero() const {
  return data == 0;
}

template<>
inline int FP_NR<long double>::is_nan() const {
  return data != data;
}

template<>
inline int FP_NR<long double>::is_finite() const {
  return isfinite(data);
}

/* arithmetic */
template<>
inline void FP_NR<long double>::add(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t /*rnd*/) {
  data = b.data + c.data;
}

template<>
inline void FP_NR<long double>::sub(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t /*rnd*/) {
  data = b.data - c.data;
}

template<>
inline void FP_NR<long double>::mul(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t /*rnd*/) {
  data = b.data * c.data;
}

template<>
inline void FP_NR<long double>::mul_d(const FP_NR<long double>& b, const double c, mp_rnd_t /*rnd*/) {
  data = b.data * c;
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
inline void FP_NR<long double>::addmul(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t /*rnd*/) {
  data = data + b.data * c.data;
}

template<>
inline void FP_NR<long double>::submul(const FP_NR<long double>& b, const FP_NR<long double>& c, mp_rnd_t /*rnd*/) {
  data = data - b.data * c.data;
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
inline void FP_NR<long double>::sqrt(const FP_NR<long double>& s, mp_rnd_t /*rnd*/) {
  data = sqrtl(s.data);
}

template<>
inline void FP_NR<long double>::root(const FP_NR<long double>& a, unsigned int k, mp_rnd_t /*rnd*/) {
  data = powl(a.data, 1/(static_cast<double>(k)));
}

template<>
inline void FP_NR<long double>::neg(const FP_NR<long double>& b) {
  data = -b.data;
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
inline void FP_NR<long double>::rnd_we(const FP_NR<long double>& b, long expo_add) {
  // If data == 0.0, exponent() is undefined, but both branches will work
  if (b.exponent() + expo_add >= numeric_limits<long double>::digits)
    data = b.data;
  else
    data = ldexpl(rintl(ldexpl(b.data, expo_add)), -expo_add);
}

template<>
inline void FP_NR<long double>::floor(const FP_NR<long double>& b) {
  data = floorl(b.data);
}

template<>
inline void FP_NR<long double>::set_nan() {
  data = NAN;
}

template<>
inline void FP_NR<long double>::swap(FP_NR<long double>& a) {
  std::swap(data, a.data);
}

template<>
inline void FP_NR<long double>::hypot(const FP_NR<long double>& a, const FP_NR<long double>& b, mp_rnd_t /*rnd*/) {
  data = hypotl(a.data, b.data);
}


#endif // End FPLLL_WITH_LONG_DOUBLE

FPLLL_END_NAMESPACE

#endif
