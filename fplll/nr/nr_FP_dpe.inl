/****************************
 *  F=dpe_t specialization
 ****************************/

#ifndef FPLLL_NR_FP_DPE_H
#define FPLLL_NR_FP_DPE_H

FPLLL_BEGIN_NAMESPACE

/* DPE specialization if defined */
#ifdef FPLLL_WITH_DPE

/* constructor */
template <> inline FP_NR<dpe_t>::FP_NR() { dpe_init(data); }

template <> inline FP_NR<dpe_t>::FP_NR(const FP_NR<dpe_t> &f)
{
  dpe_init(data);
  dpe_set(data, f.data);
}

template <> inline FP_NR<dpe_t>::~FP_NR() { dpe_clear(data); }

template <> inline unsigned int FP_NR<dpe_t>::get_prec() { return DPE_BITSIZE; }

template <> inline unsigned int FP_NR<dpe_t>::set_prec(unsigned int /*prec*/)
{
  return get_prec();  // ignored
}

/* return data */
template <> inline double FP_NR<dpe_t>::get_d(mp_rnd_t /*rnd*/) const { return dpe_get_d(data); }

template <> inline void FP_NR<dpe_t>::get_mpfr(mpfr_t r, mp_rnd_t rnd) const
{
  dpe_get_mpfr(r, data, rnd);
}

template <> inline void FP_NR<dpe_t>::set_mpfr(mpfr_t r, mp_rnd_t /*rnd*/)
{
  dpe_set_mpfr(data, r);
}

template <> inline long FP_NR<dpe_t>::get_si() const { return dpe_get_si(data); }

template <> inline long FP_NR<dpe_t>::exponent() const { return DPE_EXP(data); }

template <> inline long FP_NR<dpe_t>::get_si_exp(long &expo) const
{
  long result;
  expo = dpe_get_si_exp(&result, data);
  if (dpe_zero_p(data))
    expo = 0;
  else if (expo < 0)
  {
    /* NOTE: conversion of result to double is exact even if
        sizeof(long) = 8 */
    result = static_cast<long>(ldexp(static_cast<double>(result), expo));
    expo   = 0;
  }
  return result;
}

template <> inline long FP_NR<dpe_t>::get_si_exp_we(long &expo, long /*expo_add*/) const
{
  return get_si_exp(expo);  // NOTE: expo_add = 0
}

/*  comparison */
template <> inline int FP_NR<dpe_t>::cmp(const FP_NR<dpe_t> &a) const
{
  return dpe_cmp(data, a.data);
}

template <> inline int FP_NR<dpe_t>::cmp(double a) const { return dpe_cmp_d(data, a); }

template <> inline int FP_NR<dpe_t>::sgn() const { return cmp(0.0); }

/* operators */
template <> inline FP_NR<dpe_t> &FP_NR<dpe_t>::operator=(const FP_NR<dpe_t> &f)
{
  dpe_set(data, f.data);
  return *this;
}

template <> inline FP_NR<dpe_t> &FP_NR<dpe_t>::operator=(double d)
{
  dpe_set_d(data, d);
  return *this;
}

template <> inline FP_NR<dpe_t> &FP_NR<dpe_t>::operator=(const char *s)
{
  DPE_DOUBLE d;
#ifdef DPE_USE_DOUBLE
  d = strtod(s, NULL);
#else
  d = strtold(s, NULL);
#endif
  dpe_set_d(data, d);
  return *this;
}

template <> inline bool FP_NR<dpe_t>::operator<=(const FP_NR<dpe_t> &a) const
{
  return dpe_cmp(data, a.data) <= 0;
}

template <> inline bool FP_NR<dpe_t>::operator<=(double a) const { return dpe_cmp_d(data, a) <= 0; }

template <> inline bool FP_NR<dpe_t>::operator>=(const FP_NR<dpe_t> &a) const
{
  return dpe_cmp(data, a.data) >= 0;
}

template <> inline bool FP_NR<dpe_t>::operator>=(double a) const { return dpe_cmp_d(data, a) >= 0; }

template <> inline bool FP_NR<dpe_t>::operator<(const FP_NR<dpe_t> &a) const
{
  return dpe_cmp(data, a.data) < 0;
}

template <> inline bool FP_NR<dpe_t>::operator<(double a) const { return dpe_cmp_d(data, a) < 0; }

template <> inline bool FP_NR<dpe_t>::operator>(const FP_NR<dpe_t> &a) const
{
  return dpe_cmp(data, a.data) > 0;
}

template <> inline bool FP_NR<dpe_t>::operator>(double a) const { return dpe_cmp_d(data, a) > 0; }

template <> inline bool FP_NR<dpe_t>::is_zero() const { return dpe_zero_p(data); }

template <> inline int FP_NR<dpe_t>::is_nan() const { return DPE_MANT(data) != DPE_MANT(data); }

template <> inline int FP_NR<dpe_t>::is_finite() const { return isfinite(DPE_MANT(data)); }

/* arithmetic */
template <>
inline void FP_NR<dpe_t>::add(const FP_NR<dpe_t> &a, const FP_NR<dpe_t> &b, mp_rnd_t /*rnd*/)
{
  dpe_add(data, a.data, b.data);
}

template <>
inline void FP_NR<dpe_t>::sub(const FP_NR<dpe_t> &a, const FP_NR<dpe_t> &b, mp_rnd_t /*rnd*/)
{
  dpe_sub(data, a.data, b.data);
}

template <>
inline void FP_NR<dpe_t>::mul(const FP_NR<dpe_t> &a, const FP_NR<dpe_t> &b, mp_rnd_t /*rnd*/)
{
  dpe_mul(data, a.data, b.data);
}

template <> inline void FP_NR<dpe_t>::mul_d(const FP_NR<dpe_t> &a, const double b, mp_rnd_t /*rnd*/)
{
  dpe_mul_d(data, a.data, b);
}

template <> inline void FP_NR<dpe_t>::mul_2si(const FP_NR<dpe_t> &a, long b)
{
  dpe_mul_2si(data, a.data, b);
}

template <>
inline void FP_NR<dpe_t>::div(const FP_NR<dpe_t> &a, const FP_NR<dpe_t> &b, mp_rnd_t /*rnd*/)
{
  dpe_div(data, a.data, b.data);
}

template <>
inline void FP_NR<dpe_t>::addmul(const FP_NR<dpe_t> &a, const FP_NR<dpe_t> &b, mp_rnd_t /*rnd*/)
{
  dpe_addmul(data, a.data, b.data);
}

template <>
inline void FP_NR<dpe_t>::submul(const FP_NR<dpe_t> &a, const FP_NR<dpe_t> &b, mp_rnd_t /*rnd*/)
{
  dpe_submul(data, a.data, b.data);
}

template <> inline void FP_NR<dpe_t>::pow_si(const FP_NR<dpe_t> &a, long b, mp_rnd_t /*rnd*/)
{
  dpe_pow_si(data, a.data, b);
}

template <> inline void FP_NR<dpe_t>::exponential(const FP_NR<dpe_t> &a, mp_rnd_t /*rnd*/)
{
  // dpe_ugly_exp(data, a.data);
  dpe_exponential(data, a.data);
}

template <> inline void FP_NR<dpe_t>::log(const FP_NR<dpe_t> &a, mp_rnd_t /*rnd*/)
{
  dpe_log(data, a.data);
  // dpe_ugly_log(data, a.data);
}

template <> inline void FP_NR<dpe_t>::sqrt(const FP_NR<dpe_t> &a, mp_rnd_t /*rnd*/)
{
  dpe_sqrt(data, a.data);
}

/* root() function not available */

template <> inline void FP_NR<dpe_t>::neg(const FP_NR<dpe_t> &a) { dpe_neg(data, a.data); }

template <> inline void FP_NR<dpe_t>::abs(const FP_NR<dpe_t> &a) { dpe_abs(data, a.data); }

template <> inline void FP_NR<dpe_t>::rnd(const FP_NR<dpe_t> &a) { dpe_round(data, a.data); }

template <> inline void FP_NR<dpe_t>::rnd_we(const FP_NR<dpe_t> &a, long /*expo_add*/)
{
  dpe_round(data, a.data);  // NOTE: expo_add = 0
}

template <> inline void FP_NR<dpe_t>::floor(const FP_NR<dpe_t> &a) { dpe_floor(data, a.data); }

template <> inline void FP_NR<dpe_t>::ceil(const FP_NR<dpe_t> &a) { dpe_ceil(data, a.data); }

template <> inline void FP_NR<dpe_t>::set_nan()
{
  // dpe_set_d(data, NAN); // DPE_UNLIKELY branch in dpe_normalize
  DPE_MANT(data) = NAN;
}

template <> inline void FP_NR<dpe_t>::swap(FP_NR<dpe_t> &a) { dpe_swap(data, a.data); }

template <>
inline void FP_NR<dpe_t>::hypot(const FP_NR<dpe_t> &a, const FP_NR<dpe_t> &b, mp_rnd_t /*rnd*/)
{
  dpe_t tmp1;
  dpe_t one;
  dpe_set_d(one, 1.0);
  dpe_abs(data, a.data);  // this <- abs(a)
  dpe_abs(tmp1, b.data);  // tmp1 <- abs(b)
  if (dpe_cmp(data, tmp1))
  {                             // if this = abs_a <= abs_b = tmp1
    dpe_div(data, data, tmp1);  // this <- abs_a/abs_b = this/tmp1
    dpe_mul(data, data, data);  // this <- (abs_a/abs_b)^2
    dpe_add(data, data, one);   // this <- 1 + (abs_a/abs_b)^2
    dpe_sqrt(data, data);       // this <- sqrt(1 + (abs_a/abs_b)^2)
    dpe_mul(data, data, tmp1);  // this <- sqrt(1 + (abs_a/abs_b)^2)*abs_b
  }
  else
  {
    dpe_div(tmp1, tmp1, data);  // tmp1 <- abs_b/abs_a = tmp1/this
    dpe_mul(tmp1, tmp1, tmp1);  // tmp1 <- (abs_b/abs_a)^2
    dpe_add(tmp1, tmp1, one);   // tmp1 <- 1 + (abs_a/abs_b)^2
    dpe_sqrt(tmp1, tmp1);       // tmp1 <- sqrt(1 + (abs_b/abs_a)^2)
    dpe_mul(data, data, tmp1);  // this <- abs)a*sqrt(1 + (abs_b/abs_a)^2)
  }

  // Naive implementation:
  //
  // dpe_mul(data,a.data,a.data);
  // dpe_addmul(data,b.data,b.data);
  // dpe_sqrt(data,data);
}

/* operators FP_NR<dpe_t> */
template <> inline ostream &operator<<(ostream &os, const FP_NR<dpe_t> &x)
{
  double m = DPE_MANT(x.get_data());
  if (!isfinite(m))
    os << m;
  else
  {
    double mm = DPE_EXP(x.get_data()) * log10(2.0);
    long e10  = static_cast<long>(mm);
    m *= pow(10.0, mm - e10);
    while (m != 0 && fabs(m) < 1)
    {
      m *= 10;
      e10--;
    }
    os << m;
    if (m != 0 && e10 != 0)
      os << "e" << e10;
  }
  return os;
}

#endif  // #ifdef FPLLL_WITH_DPE

FPLLL_END_NAMESPACE

#endif
