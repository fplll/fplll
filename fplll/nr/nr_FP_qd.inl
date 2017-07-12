#ifndef FPLLL_NR_FP_QD_H
#define FPLLL_NR_FP_QD_H
/**********************************
 *  F=qd_real specialization
 *********************************/

#include "../defs.h"
#include "nr_FP.inl"
#include <qd/qd_real.h>

FPLLL_BEGIN_NAMESPACE

/* QD specialization if defined QD */


/* constructor */
template<>
inline FP_NR<qd_real>::FP_NR() {}

template<>
inline FP_NR<qd_real>::FP_NR(const FP_NR<qd_real>& f) {data = f.data;}

template<>
inline FP_NR<qd_real>::~FP_NR() {}

template<>
inline unsigned int FP_NR<qd_real>::get_prec() {
  return PREC_QD;
}

template<>
inline unsigned int FP_NR<qd_real>::set_prec(unsigned int) {
  return get_prec(); // ignored
}

/* return data */
template<>
inline double FP_NR<qd_real>::get_d(mp_rnd_t /*rnd*/) const {
  return ::to_double(data);
}

template<>
inline void FP_NR<qd_real>::get_mpfr(mpfr_t r, mp_rnd_t rnd) const {
  mpfr_set_prec (r, get_prec());
  mpfr_set_d (r, data.x[0], rnd);
  mpfr_add_d (r, r, data.x[1], rnd);
  mpfr_add_d (r, r, data.x[2], rnd);
  mpfr_add_d (r, r, data.x[3], rnd);
}

template<>
inline void FP_NR<qd_real>::set_mpfr(mpfr_t r, mp_rnd_t rnd) {
  data = mpfr_get_ld (r, rnd);
}

template<>
inline long FP_NR<qd_real>::get_si() const {
  return to_int(data);
}

template<>
inline long FP_NR<qd_real>::exponent() const {
  return ilogb(::to_double(data)) + 1;
}

template<>
inline long FP_NR<qd_real>::get_si_exp_we(long& expo, long expo_add) const {
  if (data == 0)
    expo = 0;
  else
    expo = max(exponent() + expo_add - numeric_limits<long>::digits, 0L);
  return static_cast<long>((::ldexp(data, expo_add - expo)).x[0]);
}

template<>
inline long FP_NR<qd_real>::get_si_exp(long& expo) const {
  return get_si_exp_we(expo, 0);
}

/*  comparison */
template<>
inline int  FP_NR<qd_real>::cmp(const FP_NR<qd_real>& b) const {
  if (data > b.data) return 1;
  if (data < b.data) return -1;
  return 0;
}

template<>
inline int  FP_NR<qd_real>::cmp(double b) const {
  if (data > b) return 1;
  if (data < b) return -1;
  return 0;
}

template<>
inline int FP_NR<qd_real>::sgn() const {
  if (data > 0) return 1;
  if (data < 0) return -1;
  return 0;
}

/* operators */
template<>
inline FP_NR<qd_real>& FP_NR<qd_real>::operator=(const FP_NR<qd_real>& f) {
  data = f.data;
  return *this;
}

template<>
inline FP_NR<qd_real>& FP_NR<qd_real>::operator=(double d) {
  data = d;
  return *this;
}

template<>
inline FP_NR<qd_real>& FP_NR<qd_real>::operator=(const char* s) {
  data.read(s, data);
  return *this;
}

template<>
inline bool FP_NR<qd_real>::operator<=(const FP_NR<qd_real>& a) const {
  return data <= a.data;
}

template<>
inline bool FP_NR<qd_real>::operator<=(double a) const {
  return data <= a;
}

template<>
inline bool FP_NR<qd_real>::operator>=(const FP_NR<qd_real>& a) const {
  return data >= a.data;
}

template<>
inline bool FP_NR<qd_real>::operator>=(double a) const {
  return data >= a;
}

template<>
inline bool FP_NR<qd_real>::operator<(const FP_NR<qd_real>& a) const {
  return data < a.data;
}

template<>
inline bool FP_NR<qd_real>::operator<(double a) const {
  return data < a;
}

template<>
inline bool FP_NR<qd_real>::operator>(const FP_NR<qd_real>& a) const {
  return data > a.data;
}

template<>
inline bool FP_NR<qd_real>::operator>(double a) const {
  return data > a;
}

template<>
inline bool FP_NR<qd_real>::is_zero() const {
  return data == 0;
}

template<>
inline int FP_NR<qd_real>::is_nan() const {
  return data.isnan();
}

template<>
inline int FP_NR<qd_real>::is_finite() const {
  return ::isfinite(data);
}

/* arithmetic */
template<>
inline void FP_NR<qd_real>::add(const FP_NR<qd_real>& b, const FP_NR<qd_real>& c, mp_rnd_t /*rnd*/) {
  data = b.data + c.data;
}

template<>
inline void FP_NR<qd_real>::sub(const FP_NR<qd_real>& b, const FP_NR<qd_real>& c, mp_rnd_t /*rnd*/) {
  data = b.data - c.data;
}

template<>
inline void FP_NR<qd_real>::mul(const FP_NR<qd_real>& b, const FP_NR<qd_real>& c, mp_rnd_t /*rnd*/) {
  data = b.data * c.data;
}

template<>
inline void FP_NR<qd_real>::mul_d(const FP_NR<qd_real>& b, const double c, mp_rnd_t /*rnd*/) {
  data = b.data * c;
}

template<>
inline void FP_NR<qd_real>::mul_2si(const FP_NR<qd_real>& b, long c) {
  data = ::ldexp(b.data, static_cast<int>(c));
}

template<>
inline void FP_NR<qd_real>::div(const FP_NR<qd_real>& a, const FP_NR<qd_real>& b, mp_rnd_t rnd) {
  data = a.data / b.data;
}

template<>
inline void FP_NR<qd_real>::addmul(const FP_NR<qd_real>& b, const FP_NR<qd_real>& c, mp_rnd_t /*rnd*/) {
  data = data + b.data * c.data;
}

template<>
inline void FP_NR<qd_real>::submul(const FP_NR<qd_real>& b, const FP_NR<qd_real>& c, mp_rnd_t /*rnd*/) {
  data = data - b.data * c.data;
}

template<>
inline void FP_NR<qd_real>::pow_si(const FP_NR<qd_real>& a, long b, mp_rnd_t /*rnd*/) {
  data = ::pow(a.data, static_cast<int>(b));
}

template<>
inline void FP_NR<qd_real>::exponential(const FP_NR<qd_real>& a, mp_rnd_t /*rnd*/) {
  data = ::exp(a.data);
}

template<>
inline void FP_NR<qd_real>::log(const FP_NR<qd_real>& a, mp_rnd_t /*rnd*/) {
  data = ::log(a.data);
}

template<>
inline void FP_NR<qd_real>::sqrt(const FP_NR<qd_real>& a, mp_rnd_t /*rnd*/) {
  data = ::sqrt(a.data);
}

template<>
inline void FP_NR<qd_real>::root(const FP_NR<qd_real>& a, unsigned int k, mp_rnd_t /*rnd*/) {
  data = ::nroot(a.data, k);
}

template<>
inline void FP_NR<qd_real>::neg(const FP_NR<qd_real>& b) {
  data = -b.data;
}

template<>
inline void FP_NR<qd_real>::abs(const FP_NR<qd_real>& b) {
  data = ::abs(b.data);
}

template<>
inline void FP_NR<qd_real>::rnd(const FP_NR<qd_real>& b) {
  data = ::nint(b.data);
}

template<>
inline void FP_NR<qd_real>::rnd_we(const FP_NR<qd_real>& b, long expo_add) {
  /* quad-double has same expo limit as double format*/
  if (b.exponent() + expo_add >= numeric_limits<double>::digits)
    data = b.data;
  else
    data = ::ldexp(::nint(::ldexp(b.data, expo_add)), -expo_add);
}

template<>
inline void FP_NR<qd_real>::floor(const FP_NR<qd_real>& b) {
  data = ::floor(b.data);
}

template<>
inline void FP_NR<qd_real>::set_nan() {
  data = NAN;
}

template<>
inline void FP_NR<qd_real>::swap(FP_NR<qd_real>& a) {
  qd_real t = a.data;
  a.data = data;
  data = t;
}

template<>
inline void FP_NR<dd_real>::hypot(const FP_NR<qd_real>& a, const FP_NR<qd_real>& b, mp_rnd_t rnd) {
  // Maybe decrease temporary 
  // variables with one, by
  // putting the absolute value of a
  // into this.data

    dd_real abs_a = ::abs(a.data);
    dd_real abs_b = ::abs(b.data);
    
    if (abs_a > abs_b) {
        data = abs_a* ::sqrt(1.0 + (::pow(abs_b/abs_a,2)));
    } else {
        data =  abs_b* ::sqrt(1.0 + (::pow(abs_a/abs_b,2)));   
    }
    
}

FPLLL_END_NAMESPACE

#endif
