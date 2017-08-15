/****************************
 *  F=mpfr_t specialization
 ****************************/

#ifndef FPLLL_NR_FP_MPFR_H
#define FPLLL_NR_FP_MPFR_H

FPLLL_BEGIN_NAMESPACE

/* MPFR specialization */

/* constructor */
template<>
inline FP_NR<>::FP_NR() {mpfr_init(data);}

template<>
inline FP_NR<>::FP_NR(const FP_NR<>& f) {
  mpfr_init_set(data, f.data, GMP_RNDN);
}

template<>
inline FP_NR<>::~FP_NR() {mpfr_clear(data);}

template<>
inline unsigned int FP_NR<>::get_prec() {
  return mpfr_get_default_prec();
}

template<>
inline unsigned int FP_NR<>::set_prec(unsigned int prec) {
  int old_prec = get_prec();
  mpfr_set_default_prec(prec);
  return old_prec;
}

/* return data */
template<>
inline double FP_NR<>::get_d(mp_rnd_t rnd) const {
  return mpfr_get_d(data, rnd);
}

template<>
inline void FP_NR<>::get_mpfr(mpfr_t r, mp_rnd_t rnd) const {
  mpfr_set(r, data, rnd);
}

template<>
inline void FP_NR<>::set_mpfr(mpfr_t r, mp_rnd_t rnd) {
  mpfr_set(data, r, rnd);
}

template<>
inline long FP_NR<>::get_si() const {
  return mpfr_get_si(data,GMP_RNDN);
}

template<>
inline long FP_NR<>::exponent() const {
  return mpfr_get_exp(data);
}

template<>
inline long FP_NR<>::get_si_exp(long& expo) const {
  if (mpfr_zero_p(data)) {
    expo = 0;
  }
  else {
    expo = max(exponent() - numeric_limits<long>::digits, 0L);
  }
  mpfr_t& nc_data = const_cast<mpfr_t&>(data);
  mpfr_div_2si(nc_data, nc_data, expo, GMP_RNDN);
  long result = mpfr_get_si(nc_data, GMP_RNDZ);
  mpfr_mul_2si(nc_data, nc_data, expo, GMP_RNDN);
  return result;
}

template<>
inline long FP_NR<>::get_si_exp_we(long& expo, long /*expo_add*/) const {
  return get_si_exp(expo);  // NOTE: expo_add = 0
}

/*  comparison */
template<>
inline int FP_NR<>::cmp(const FP_NR<>& a) const {
  return mpfr_cmp(data, a.data);
}

template<>
inline int FP_NR<>::cmp(double a) const {
  return mpfr_cmp_d(data, a);
}

template<>
inline int FP_NR<>::sgn() const {
  return mpfr_sgn(data);
}

/*operators */
template<>
inline FP_NR<>& FP_NR<>::operator=(const FP_NR<>& a) {
  mpfr_set(data, a.data, GMP_RNDN);
  return *this;
}

template<>
inline FP_NR<>& FP_NR<>::operator=(double a) {
  mpfr_set_d(data, a, GMP_RNDN);
  return *this;
}

template<>
inline FP_NR<>& FP_NR<>::operator=(const char *s) {
  mpfr_set_str(data, s, 10, GMP_RNDN);
  return *this;
}

template<>
inline FP_NR<>& FP_NR<>::operator=(mpfr_t& a) {
  mpfr_set(data, a, GMP_RNDN);
  return *this;
}

template<>
inline bool FP_NR<>::operator<=(const FP_NR<>& a) const {
  return mpfr_cmp(data, a.data) <= 0;
}

template<>
inline bool FP_NR<>::operator<=(double a) const {
  return mpfr_cmp_d(data, a) <= 0;
}

template<>
inline bool FP_NR<>::operator>=(const FP_NR<>& a) const {
  return mpfr_cmp(data, a.data) >= 0;
}

template<>
inline bool FP_NR<>::operator>=(double a) const {
  return mpfr_cmp_d(data, a) >= 0;
}

template<>
inline bool FP_NR<>::operator<(const FP_NR<>& a) const {
  return mpfr_cmp(data, a.data) < 0;
}

template<>
inline bool FP_NR<>::operator<(double a) const {
  return mpfr_cmp_d(data, a) < 0;
}

template<>
inline bool FP_NR<>::operator>(const FP_NR<>& a) const {
  return mpfr_cmp(data, a.data) > 0;
}

template<>
inline bool FP_NR<>::operator>(double a) const {
  return mpfr_cmp_d(data, a) > 0;
}

template<>
inline bool FP_NR<>::is_zero() const {
  return mpfr_zero_p(data);
}

template<>
inline int FP_NR<>::is_nan() const {
  return mpfr_nan_p(data);
}

template<>
inline int FP_NR<>::is_finite() const {
  return mpfr_number_p(data);
}

/* arithmetic */
template<>
inline void FP_NR<>::add(const FP_NR<>& a, const FP_NR<>& b, mp_rnd_t rnd) {
  mpfr_add(data, a.data, b.data, rnd);
}

template<>
inline void FP_NR<>::sub(const FP_NR<>& a, const FP_NR<>& b, mp_rnd_t rnd) {
  mpfr_sub(data, a.data, b.data, rnd);
}

template<>
inline void FP_NR<>::mul(const FP_NR<>& a, const FP_NR<>& b, mp_rnd_t rnd) {
  mpfr_mul(data, a.data, b.data, rnd);
}

template<>
inline void FP_NR<>::mul_d(const FP_NR<>& a, double b, mp_rnd_t rnd) {
  mpfr_mul_d(data, a.data, b, rnd);
}

template<>
inline void FP_NR<>::mul_mpfr(const FP_NR<>& a, const mpfr_t b, mp_rnd_t rnd) {
  mpfr_mul(data, a.data, b, rnd);
}

template<>
inline void FP_NR<>::mul_2si(const FP_NR<>& a, long b) {
  mpfr_mul_2si(data, a.data, b, GMP_RNDN);
}

template<>
inline void FP_NR<>::div(const FP_NR<>& a, const FP_NR<>& b, mp_rnd_t rnd) {
  mpfr_div(data, a.data, b.data, rnd);
}

template<>
inline void FP_NR<>::div_d(const FP_NR<>& a, const double b, mp_rnd_t rnd) {
  mpfr_div_d(data, a.data, b, rnd);
}

template<>
inline void FP_NR<>::addmul(const FP_NR<>& b, const FP_NR<>& c, mp_rnd_t rnd) {
  mpfr_fma(data, b.data, c.data, data, rnd);
}

template<>
inline void FP_NR<>::submul(const FP_NR<>& b, const FP_NR<>& c, mp_rnd_t rnd) {
  mpfr_fms(data, b.data, c.data, data, rnd);
  mpfr_neg(data, data, GMP_RNDN); // Exact
}

template<>
inline void FP_NR<>::pow_si(const FP_NR<>& a, long b, mp_rnd_t rnd) {
  mpfr_pow_si(data, a.data, b, rnd);
}

template<>
inline void FP_NR<>::exponential(const FP_NR<>& a, mp_rnd_t rnd) {
  mpfr_exp(data, a.data, rnd);
}

template<>
inline void FP_NR<>::log(const FP_NR<>& a, mp_rnd_t rnd) {
  mpfr_log(data, a.data, rnd);
}

template<>
inline void FP_NR<>::sqrt(const FP_NR<>& a, mp_rnd_t rnd) {
  mpfr_sqrt(data, a.data, rnd);
}

template<>
inline void FP_NR<>::root(const FP_NR<>& a, unsigned int k, mp_rnd_t rnd) {
  mpfr_root(data, a.data, k, rnd);
}

template<>
inline void FP_NR<>::neg(const FP_NR<>& a) {
  mpfr_neg(data, a.data, GMP_RNDN);
}

template<>
inline void FP_NR<>::abs(const FP_NR<>& a) {
  mpfr_abs(data, a.data, GMP_RNDN);
}
template<>
inline void FP_NR<>::rnd(const FP_NR<>& a) {
  mpfr_round(data, a.data);
}
template<>
inline void FP_NR<>::rnd_we(const FP_NR<>& a, long /*expo_add*/) {
  mpfr_round(data, a.data); // NOTE: expo_add = 0
}
template<>
inline void FP_NR<>::floor(const FP_NR<>& a) {
  mpfr_floor(data, a.data);
}

template<>
inline void FP_NR<>::set_nan() {
  mpfr_set_nan(data);
}

template<>
inline void FP_NR<>::swap(FP_NR<>& a) {
  mpfr_swap(data, a.data);
}

/* the hypot function for Givens rotations */
template<>
inline void FP_NR<mpfr_t>::hypot(const FP_NR<mpfr_t>& a, const FP_NR<mpfr_t>& b, mp_rnd_t rnd) {
  mpfr_hypot(data, a.data, b.data, rnd);
}




/* operators FP_NR<> */
template<>
inline ostream& operator<<(ostream& os, const FP_NR<>& x) {
  mp_exp_t e;
  char* s = mpfr_get_str(NULL, &e, 10, os.precision(), x.get_data(), GMP_RNDN);
  char* p = s;
  if (*p == '-') {
    os << *p;
    p++;
  }
  if (*p == '@' || *p == 0)
    os << p;
  else if (*p == '0')
    os << *p;
  else {
    os << *p << '.' << p + 1;
    if (e - 1 != 0) os << 'e' << e - 1;
  }
  mpfr_free_str(s);
  return os;
}


FPLLL_END_NAMESPACE

#endif
