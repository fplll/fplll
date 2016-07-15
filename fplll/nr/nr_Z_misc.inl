/************************
 * Extra utils for NR_Z
 ************************/

#ifndef FPLLL_NR_Z_MISC_H
#define FPLLL_NR_Z_MISC_H

FPLLL_BEGIN_NAMESPACE


/* -----------------------------
 *   get_f_exp (Z_NR --> FP_NR)
 * ----------------------------- */

/* get_f_exp (long to double) */
#ifdef FPLLL_WITH_ZLONG

template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<double>& f, long& expo) {
  int intExpo;
  f.get_data() = frexp(static_cast<double>(data), &intExpo);
  expo = intExpo;
}

#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<long double>& f, long& expo) {
  int intExpo;
  f.get_data() = frexpl(static_cast<long double>(data), &intExpo);
  expo = intExpo;
}
#endif

#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<dpe_t>& f, long& expo) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for dpe_t");
}
#endif

#ifdef FPLLL_WITH_QD
template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<dd_real>& f, long& expo) {
  int intExpo;
  f.get_data() = frexp(static_cast<double>(data), &intExpo);
  expo = intExpo;
}
template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<qd_real>& f, long& expo) {
  int intExpo;
  f.get_data() = frexp(static_cast<double>(data), &intExpo);
  expo = intExpo;
}
#endif

template<> template<>
inline void Z_NR<long>::get_f_exp(FP_NR<mpfr_t>& /*f*/, long& /*expo*/) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for mpfr_t");
}

#endif 


/* get_f_exp (double to double) */
#ifdef FPLLL_WITH_ZDOUBLE

template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<double>& f, long& expo) {
  int intExpo;
  f.get_data() = frexp(data, &intExpo);
  expo = intExpo;
}

#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<long double>& f, long& expo) {
  int intExpo;
  f.get_data() = static_cast<long double>(frexp(data, &intExpo));
  expo = intExpo;
}
#endif

#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<dpe_t>& /*f*/, long& /*expo*/) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for dpe_t");
}
#endif

#ifdef FPLLL_WITH_QD
template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<dd_real>& f, long& expo) {
  int intExpo;
  f.get_data() = frexp(data, &intExpo);
  expo = intExpo;
}
template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<qd_real>& f, long& expo) {
  int intExpo;
  f.get_data() = frexp(data, &intExpo);
  expo = intExpo;
}
#endif

template<> template<>
inline void Z_NR<double>::get_f_exp(FP_NR<mpfr_t>& f, long& expo) {
  int intExpo;
  f = frexp(data, &intExpo);
  expo = intExpo;
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for mpfr_t");
}

#endif


/* get_f_exp (mpz_t to double) */
template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<double>& f, long& expo) {
  f.get_data() = mpz_get_d_2exp(&expo, data);
}

#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<long double>& f, long& expo) {
  f.get_data() = LDConvHelper::mpz_get_ld_2exp(&expo, data);
}
#endif

#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<dpe_t>& /*f*/, long& /*expo*/) {
  FPLLL_DEBUG_ABORT("get_f_exp unimplemented for dpe_t");
}
#endif

#ifdef FPLLL_WITH_QD
template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<dd_real>& f, long& expo) {
  f.get_data() = mpz_get_d_2exp(&expo, data);
}
template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<qd_real>& f, long& expo) {
  f.get_data() = mpz_get_d_2exp(&expo, data);
}
#endif

template<> template<>
inline void Z_NR<mpz_t>::get_f_exp(FP_NR<mpfr_t>& f, long& expo) {
  f = mpz_get_d_2exp(&expo, data);
}


/* -------------------------
 *   set_f (FP_NR --> Z_NR)
 * ------------------------- */

/* set_f (FP_NR to long) */
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

#ifdef FPLLL_WITH_QD
template<> template<>
inline void Z_NR<long>::set_f(const FP_NR<dd_real>& a) {
  data = a.get_si();
}

template<> template<>
inline void Z_NR<long>::set_f(const FP_NR<qd_real>& a) {
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

/* set_f (FP_NR to Zdouble) */
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

#ifdef FPLLL_WITH_QD
template<> template<>
inline void Z_NR<double>::set_f(const FP_NR<dd_real>& a) {
  data = a.get_d();
}

template<> template<>
inline void Z_NR<double>::set_f(const FP_NR<qd_real>& a) {
  data = a.get_d();
}
#endif

template<> template<>
inline void Z_NR<double>::set_f(const FP_NR<mpfr_t>& a) {
  data = a.get_d();
}
#endif // #ifdef FPLLL_WITH_ZDOUBLE

/* set_f (FP_NR to mpz_t) */
template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<double>& a) {
  mpz_set_d(data, a.get_data());
}

#ifdef FPLLL_WITH_LONG_DOUBLE
template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<long double>& a) {
  LDConvHelper::mpz_set_ld(data, a.get_data());
}
#endif

#ifdef FPLLL_WITH_DPE
template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<dpe_t>& a) {
  dpe_get_z(data, const_cast<dpe_t&>(a.get_data()));
}
#endif

#ifdef FPLLL_WITH_QD
template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<dd_real>& a) {
  mpfr_t t;
  mpfr_init(t);
  a.get_mpfr (t, GMP_RNDN);
  mpfr_get_z(data, t, GMP_RNDN);  
  mpfr_clear(t);
}

template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<qd_real>& a) {
  mpfr_t t;
  mpfr_init(t);
  a.get_mpfr (t, GMP_RNDN);
  mpfr_get_z(data, t, GMP_RNDN);  
  mpfr_clear(t);
}
#endif

template<> template<>
inline void Z_NR<mpz_t>::set_f(const FP_NR<mpfr_t>& a) {
  mpfr_get_z(data, a.get_data(), GMP_RNDN);
}

FPLLL_END_NAMESPACE

#endif
