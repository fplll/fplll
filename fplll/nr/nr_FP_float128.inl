/****************************
 *  F=__float128 specialization
 ****************************/

#ifndef FPLLL_NR_FP_FLOAT128_H
#define FPLLL_NR_FP_FLOAT128_H

extern "C" {
#include <quadmath.h>
}
#include "../defs.h"
#include "nr_FP.inl"

FPLLL_BEGIN_NAMESPACE


/* __float128 specialization */


/* constructor */
template<>
inline FP_NR<__float128>::FP_NR(): data(0.0q) {}

template<> 
inline FP_NR<__float128>::FP_NR(const FP_NR<__float128>& f) {data = f.data;}

template<>
inline FP_NR<__float128>::~FP_NR() {}

template<>
inline unsigned int FP_NR<__float128>::get_prec() {return PREC_FLOAT128;}

template<>
inline unsigned int FP_NR<__float128>::set_prec(unsigned int) {
  return get_prec(); // ignored
}

/* return data */
template<>
inline double FP_NR<__float128>::get_d(mp_rnd_t) const {
  return static_cast<__float128>(data);
}

template<>
inline void FP_NR<__float128>::get_mpfr(mpfr_t r, mp_rnd_t rnd) const {
  //mpfr_set_float128_local (r, data, rnd);
}

template<>
inline void FP_NR<__float128>::set_mpfr(mpfr_t r, mp_rnd_t rnd) {
  //data = mpfr_get_float128_local (r, rnd);
}

template<>
inline long FP_NR<__float128>::get_si() const {
  return truncq(data);
}

template<>
inline long FP_NR<__float128>::exponent() const {
  return ilogbq(data) + 1;
}

template<>
inline long FP_NR<__float128>::get_si_exp_we(long& expo, long expo_add) const {
  if (data == 0)
    expo = 0;
  else
    expo = max(exponent() + expo_add - numeric_limits<long>::digits, 0L);
  return static_cast<long>(ldexpq(data, expo_add - expo));
}

template<>
inline long FP_NR<__float128>::get_si_exp(long& expo) const {
  return get_si_exp_we(expo, 0);
}

/*  comparison */
template<>
inline int FP_NR<__float128>::cmp(const FP_NR<__float128>& b) const {
  if (data > b.data) return 1;
  if (data < b.data) return -1;
  return 0;
}

template<>
inline int FP_NR<__float128>::cmp(double b) const {
  if (data > b) return 1;
  if (data < b) return -1;
  return 0;
}

template<>
inline int FP_NR<__float128>::sgn() const {
  if (data > 0) return 1;
  if (data < 0) return -1;
  return 0;
}

/* operators */
template<>
inline FP_NR<__float128>& FP_NR<__float128>::operator=(const FP_NR<__float128>& f) {
  data = f.data;
  return *this;
}

template<>
inline FP_NR<__float128>& FP_NR<__float128>::operator=(double d) {
  data = d; // check this????
  return *this;
}

template<>
inline FP_NR<__float128>& FP_NR<__float128>::operator=(const char *s) {
  data = strtoflt128(s, NULL);
  return *this;
}

template<>
inline bool FP_NR<__float128>::operator<=(const FP_NR<__float128>& a) const {
  return data <= a.data;
}

template<>
inline bool FP_NR<__float128>::operator<=(double a) const {
  return data <= a;
}

template<>
inline bool FP_NR<__float128>::operator>=(const FP_NR<__float128>& a) const {
  return data >= a.data;
}

template<>
inline bool FP_NR<__float128>::operator>=(double a) const {
  return data >= a;
}

template<>
inline bool FP_NR<__float128>::operator<(const FP_NR<__float128>& a) const {
  return data < a.data;
}

template<>
inline bool FP_NR<__float128>::operator<(double a) const {
  return data < a;
}

template<>
inline bool FP_NR<__float128>::operator>(const FP_NR<__float128>& a) const {
  return data > a.data;
}

template<>
inline bool FP_NR<__float128>::operator>(double a) const {
  return data > a; // check?????
}

template<>
inline bool FP_NR<__float128>::is_zero() const {
  return data == 0.0q;
}

template<>
inline int FP_NR<__float128>::is_nan() const {
  return data != data;
}

template<>
inline int FP_NR<__float128>::is_finite() const {
  return finiteq(data);
}

/* arithmetic */
template<>
inline void FP_NR<__float128>::add(const FP_NR<__float128>& b, const FP_NR<__float128>& c, mp_rnd_t /*rnd*/) {
  data = b.data + c.data;
}

template<>
inline void FP_NR<__float128>::sub(const FP_NR<__float128>& b, const FP_NR<__float128>& c, mp_rnd_t /*rnd*/) {
  data = b.data - c.data;
}

template<>
inline void FP_NR<__float128>::mul(const FP_NR<__float128>& b, const FP_NR<__float128>& c, mp_rnd_t /*rnd*/) {
  data = b.data * c.data;
}

template<>
inline void FP_NR<__float128>::mul_d(const FP_NR<__float128>& b, const double c, mp_rnd_t /*rnd*/) {
  data = b.data * c;
}

template<>
inline void FP_NR<__float128>::mul_2si(const FP_NR<__float128>& b, long c) {
  data = ldexpq(b.data, static_cast<int>(c));
}

template<>
inline void FP_NR<__float128>::div(const FP_NR<__float128>& b, const FP_NR<__float128>& c, mp_rnd_t /*rnd*/) {
  data = b.data / c.data;
}

template<>
inline void FP_NR<__float128>::addmul(const FP_NR<__float128>& b, const FP_NR<__float128>& c, mp_rnd_t /*rnd*/) {
  data = data + b.data * c.data;
}

template<>
inline void FP_NR<__float128>::submul(const FP_NR<__float128>& b, const FP_NR<__float128>& c, mp_rnd_t /*rnd*/) {
  data = data - b.data * c.data;
}

template<>
inline void FP_NR<__float128>::pow_si(const FP_NR<__float128>& a, long b, mp_rnd_t /*rnd*/) {
  data = ::powq(a.data, static_cast<__float128>(b));
}

template<>
inline void FP_NR<__float128>::exponential(const FP_NR<__float128>& a, mp_rnd_t /*rnd*/) {
  data = ::expq(a.data);
}

template<>
inline void FP_NR<__float128>::log(const FP_NR<__float128>& a, mp_rnd_t /*rnd*/) {
  data = ::logq(a.data);
}

template<>
inline void FP_NR<__float128>::sqrt(const FP_NR<__float128>& s, mp_rnd_t /*rnd*/) {
  data = ::sqrtq(s.data);
}

template<>
inline void FP_NR<__float128>::root(const FP_NR<__float128>& a, unsigned int k, mp_rnd_t /*rnd*/) {
  data = ::powq(a.data, 1/(static_cast<__float128>(k)));
}

template<>
inline void FP_NR<__float128>::neg(const FP_NR<__float128>& b) {
  data = -b.data;
}

template<>
inline void FP_NR<__float128>::abs(const FP_NR<__float128>& b) {
  data = b.data;
  if (data < 0) data = -data;
}

template<>
inline void FP_NR<__float128>::rnd(const FP_NR<__float128>& b) {
  data = rintq(b.data);
}

template<>
inline void FP_NR<__float128>::rnd_we(const FP_NR<__float128>& b, long expo_add) {
  // If data == 0.0, exponent() is undefined, but both branches will
  // CHECK this??
  if (b.exponent() + expo_add >= 128) {
    data = b.data;
  }
  else
    data = ldexpq(rintq(ldexpq(b.data, expo_add)), -expo_add);
}

template<>
inline void FP_NR<__float128>::floor(const FP_NR<__float128>& b) {
  data = floorq(b.data);
}

template<>
inline void FP_NR<__float128>::set_nan() {
  data = NAN;
}

template<>
inline void FP_NR<__float128>::swap(FP_NR<__float128>& a) {
  std::swap(data, a.data);
}

/* hypot function for Givens */
template<>
inline void FP_NR<__float128>::hypot(const FP_NR<__float128>& a, const FP_NR<__float128>& b, mp_rnd_t /*rnd*/) {
  data = hypotq(a.data,b.data);
}


/* operators FP_NR<__float_128> */
template<>
inline ostream& operator<<(ostream& os, const FP_NR<__float128>& x) {
  char buf[128];
  quadmath_snprintf (buf, sizeof buf, "%+-#*.20Qe", 10, x.get_data());
  os << buf;
  return os;
}


FPLLL_END_NAMESPACE

#endif
