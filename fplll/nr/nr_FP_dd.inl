#ifndef FPLLL_NR_FP_DD_H
#define FPLLL_NR_FP_DD_H
/**********************************
 *  F=dd_real specialization
 *********************************/

#include "../defs.h"
#include "nr_FP.inl"
#include <qd/dd_real.h>

FPLLL_BEGIN_NAMESPACE

/* DD specialization if defined QD */

/* constructor */
template<>
inline FP_NR<dd_real>::FP_NR() {}

template<>
inline FP_NR<dd_real>::FP_NR(const FP_NR<dd_real>& f) {data = f.data;}

template<>
inline FP_NR<dd_real>::~FP_NR() {}

template<>
inline unsigned int FP_NR<dd_real>::getprec() {return PREC_DD;}

template<>
inline unsigned int FP_NR<dd_real>::setprec(unsigned int) {
  return getprec(); // ignored
}

/* return data */
template<>
inline double FP_NR<dd_real>::get_d(mp_rnd_t) const {
  return ::to_double(data);
}

template<>
inline void FP_NR<dd_real>::get_mpfr(mpfr_t r, mp_rnd_t rnd) const {
  mpfr_set_prec (r, getprec());
  mpfr_set_d (r, data._lo(), rnd);
  mpfr_add_d (r, r, data._hi(), rnd);
}

template<>
inline void FP_NR<dd_real>::set_mpfr(mpfr_t r, mp_rnd_t rnd) {
  data = mpfr_get_ld (r, rnd);
}

template<>
inline long FP_NR<dd_real>::get_si() const {
  return to_int(data);
}

template<>
inline long FP_NR<dd_real>::exponent() const {
  return ilogb(::to_double(data)) + 1;
}

template<>
inline long FP_NR<dd_real>::get_si_exp_we(long& expo, long expo_add) const {
  if (data == 0)
    expo = 0;
  else
    expo = max(exponent() + expo_add - numeric_limits<long>::digits, 0L);
  return static_cast<long>((::ldexp(data, expo_add - expo)).x[0]);
}

template<>
inline long FP_NR<dd_real>::get_si_exp(long& expo) const {
  return get_si_exp_we(expo, 0);
}

/*  comparison */
template<>
inline int  FP_NR<dd_real>::cmp(const FP_NR<dd_real>& b) const {
  if (data > b.data) return 1;
  if (data < b.data) return -1;
  return 0;
}

template<>
inline int  FP_NR<dd_real>::cmp(double b) const {
  if (data > b) return 1;
  if (data < b) return -1;
  return 0;
}

template<>
inline int FP_NR<dd_real>::sgn() const {
  if (data > 0) return 1;
  if (data < 0) return -1;
  return 0;
}

/* operators */
template<>
inline void FP_NR<dd_real>::operator=(const FP_NR<dd_real>& a) {
  data = a.data;
}

template<>
inline void FP_NR<dd_real>::operator=(double a) {
  data = a;
}

template<>
inline void FP_NR<dd_real>::operator=(const char* s) {
  data.read(s, data);
}

template<>
inline bool FP_NR<dd_real>::operator<=(const FP_NR<dd_real>& a) const {
  return data <= a.data;
}

template<>
inline bool FP_NR<dd_real>::operator<=(double a) const {
  return data <= a;
}

template<>
inline bool FP_NR<dd_real>::operator>=(const FP_NR<dd_real>& a) const {
  return data >= a.data;
}

template<>
inline bool FP_NR<dd_real>::operator>=(double a) const {
  return data >= a;
}

template<>
inline bool FP_NR<dd_real>::operator<(const FP_NR<dd_real>& a) const {
  return data < a.data;
}

template<>
inline bool FP_NR<dd_real>::operator<(double a) const {
  return data < a;
}

template<>
inline bool FP_NR<dd_real>::operator>(const FP_NR<dd_real>& a) const {
  return data > a.data;
}

template<>
inline bool FP_NR<dd_real>::operator>(double a) const {
  return data > a;
}

template<>
inline bool FP_NR<dd_real>::is_zero() const {
  return data == 0;
}

template<>
inline int FP_NR<dd_real>::is_nan() const {
  return data.isnan();
}

template<>
inline int FP_NR<dd_real>::is_finite() const {
  return data.isfinite();
}

/* arithmetic */
template<>
inline void FP_NR<dd_real>::add(const FP_NR<dd_real>& b, const FP_NR<dd_real>& c, mp_rnd_t /*rnd*/) {
  data = b.data + c.data;
}

template<>
inline void FP_NR<dd_real>::sub(const FP_NR<dd_real>& b, const FP_NR<dd_real>& c, mp_rnd_t /*rnd*/) {
  data = b.data - c.data;
}

template<>
inline void FP_NR<dd_real>::mul(const FP_NR<dd_real>& b, const FP_NR<dd_real>& c, mp_rnd_t /*rnd*/) {
  data = b.data * c.data;
}

template<>
inline void FP_NR<dd_real>::mul_d(const FP_NR<dd_real>& b, const double c, mp_rnd_t /*rnd*/) {
  data = b.data * c;
}

template<>
inline void FP_NR<dd_real>::mul_2si(const FP_NR<dd_real>& b, long c) {
  data = ::ldexp(b.data, static_cast<int>(c));
}

template<>
inline void FP_NR<dd_real>::div(const FP_NR<dd_real>& a, const FP_NR<dd_real>& b, mp_rnd_t rnd) {
  data = a.data / b.data;
}

template<>
inline void FP_NR<dd_real>::addmul(const FP_NR<dd_real>& b, const FP_NR<dd_real>& c, mp_rnd_t /*rnd*/) {
  data = data + b.data * c.data;
}

template<>
inline void FP_NR<dd_real>::submul(const FP_NR<dd_real>& b, const FP_NR<dd_real>& c, mp_rnd_t /*rnd*/) {
  data = data - b.data * c.data;
}

template<>
inline void FP_NR<dd_real>::pow_si(const FP_NR<dd_real>& a, long b, mp_rnd_t /*rnd*/) {
  data = ::pow(a.data, static_cast<int>(b));
}

template<>
inline void FP_NR<dd_real>::exponential(const FP_NR<dd_real>& a, mp_rnd_t /*rnd*/) {
  data = ::exp(a.data);
}

template<>
inline void FP_NR<dd_real>::log(const FP_NR<dd_real>& a, mp_rnd_t /*rnd*/) {
  data = ::log(a.data);
}

template<>
inline void FP_NR<dd_real>::sqrt(const FP_NR<dd_real>& a, mp_rnd_t /*rnd*/) {
  data = ::sqrt(a.data);
}

template<>
inline void FP_NR<dd_real>::root(const FP_NR<dd_real>& a, unsigned int k, mp_rnd_t /*rnd*/) {
  data = ::nroot(a.data, k);
}

template<>
inline void FP_NR<dd_real>::neg(const FP_NR<dd_real>& b) {
  data = -b.data;
}

template<>
inline void FP_NR<dd_real>::abs(const FP_NR<dd_real>& b) {
  data = ::abs(b.data);
}

template<>
inline void FP_NR<dd_real>::rnd(const FP_NR<dd_real>& b) {
  data = ::nint(b.data);
}

template<>
inline void FP_NR<dd_real>::rnd_we(const FP_NR<dd_real>& b, long expo_add) {
  /* double-double has same expo limit as double format*/
  if (b.exponent() + expo_add >= numeric_limits<double>::digits)
    data = b.data;
  else
    data = ::ldexp(::nint(::ldexp(b.data, expo_add)), -expo_add);
}

template<>
inline void FP_NR<dd_real>::floor(const FP_NR<dd_real>& b) {
  data = ::floor(b.data);
}

template<>
inline void FP_NR<dd_real>::set_nan() {
  data = NAN;
}

template<>
inline void FP_NR<dd_real>::swap(FP_NR<dd_real>& a) {
  dd_real t = a.data;
  a.data = data;
  data = t;
}


/* #ifdef FPLLL_V3_COMPAT */
#ifdef FPLLL_V3_COMPAT

template<>
inline void FP_NR<dd_real>::print() const {
  cout << data;
}

template<>
inline void FP_NR<dd_real>::printerr() const {
  cerr << data;
}

template<>
inline void FP_NR<dd_real>::set(unsigned int s) {
  data = static_cast<double>(s);
}
#endif

FPLLL_END_NAMESPACE

#endif
