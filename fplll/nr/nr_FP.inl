/*********************
 *   Float Class
 *********************/

#include "../defs.h"

#ifndef FPLLL_NR_FP_H
#define FPLLL_NR_FP_H

FPLLL_BEGIN_NAMESPACE

/* declaration Z_NR */
template <class F> class Z_NR;

/**
 * FP_NR stores floating-point numbers. This template provides a uniform
 * interface for doing floating-point computations with several underlying
 * types (double, dpe_t and mpfr_t). For all functions, the rounding
 * mode rnd is ignored unless F=mpfr_t.
 */
template <class F = mpfr_t> class FP_NR
{

  F data;

public:
  /**
   * Constructors.
   */
  inline FP_NR<F>();
  inline FP_NR<F>(const FP_NR<F> &f);
  inline ~FP_NR<F>();
  inline FP_NR<F>(const double d) : FP_NR<F>() {*this = d;};
  inline FP_NR<F>(const char *s) : FP_NR<F>() {*this = s;};

  /**
   * Returns the current precision for new FP_NR&lt;F&gt; objects.
   */
  static inline unsigned int get_prec();

  /**
   * Sets the precision of new FP_NR&lt;F&gt; objects. Returns the
     previous value. This has no effect is F != mpfr_t.
   */
  static inline unsigned int set_prec(unsigned int prec);

  /** get data */

  /**
   * Returns the internal representation of the data.
   */
  inline F &get_data() { return data; }

  /**
   * Returns the internal representation of the data.
   */
  inline const F &get_data() const { return data; }

  /**
   * Converts this object to a double. If it does not fit in a double,
   * the result is undefined.
   */
  inline double get_d(mp_rnd_t rnd = GMP_RNDN) const;

  /**
   * Convert this object to a mpfr_t r.
   */
  inline void get_mpfr(mpfr_t r, mp_rnd_t rnd = GMP_RNDN) const;

  /**
   * Convert this object from a mpfr_t r.
   */
  inline void set_mpfr(mpfr_t r, mp_rnd_t rnd = GMP_RNDN);

  /**
   * Converts this object to a long. The rounding direction is undefined.
   * If it does not fit in a long, the result is undefined.
   */
  inline long get_si() const;

  /**
   * Returns expo such value 2^(expo-1) <= value < 2^expo
   *  (expo = floor(abs(log2(abs(value)))) + 1).
   * The return value is undefined if *this == 0.
   */
  inline long exponent() const;

  /**
   * Returns x and defines expo such that
   *  trunc(value * 2^expo_add) ~= x * 2^expo
   * The '~=' is an equality if trunc(value * 2^expo_add) <= LONG_MAX.
   * expo is the minimum non-negative value such that x <= LONG_MAX.
   * expo_add must be 0 if T=dpe_t or T=mpfr_t.
   */
  inline long get_si_exp_we(long &expo, long expo_add) const;

  /**
   * Returns x and defines expo such that trunc(value) ~= x * 2^expo.
   * The '~=' is an equality if trunc(value) <= LONG_MAX.
   * expo is the minimum non-negative value such that |x| <= LONG_MAX.
   * (x is the largest among all possible |x| <= LONG_MAX).
   */
  inline long get_si_exp(long &expo) const;

  /**
   * Computes a and expo such that trunc(value) ~= a * 2^expo.
   *  The '~=' is an equality if Z=mpz_t. expo is always non-negative.
   *  Note that expo might be very close to LONG_MAX when value = 0.
   *  (in nr_FP_misc.h)
   */
  template <class Z> inline void get_z_exp(Z_NR<Z> &a, long &expo) const;

  /**
   * Computes a and expo such that trunc(value) * 2^expo_add ~= a * 2^expo.
   * The '~=' is an equality if Z=mpz_t. expo is always non-negative.
   * expo_add must be 0 if T=dpe_t or T=mpfr_t.
   *  (in nr_FP_misc.h)
   */
  template <class Z> inline void get_z_exp_we(Z_NR<Z> &a, long &expo, long expo_add) const;

  /** set data */

  /**
   * Sets the value to z.
   * (in nr_FP_misc.h)
   */
  template <class Z> inline void set_z(const Z_NR<Z> &z, mp_rnd_t rnd = GMP_RNDN);

  /** comparison */

  /**
   * 3-way comparison. Returns a positive number if *this > b,
   * a negative number if *this < b or zero is *this == b.
   */
  inline int cmp(const FP_NR<F> &b) const;

  /**
   * 3-way comparison. Returns a positive number if *this > d,
   * a negative number if *this < d or zero is *this == d.
   */
  inline int cmp(double d) const;

  /**
   * Sign. Returns a positive number, a negative number or zero if the
   * value of this object is respectively positive, negative or null.
   */
  inline int sgn() const;

  /**
   * Operator
   */
  inline FP_NR<F>& operator=(const FP_NR<F> &a);
  inline FP_NR<F>& operator=(const char *s);
  inline FP_NR<F>& operator=(double a);
  // inline FP_NR<F>& operator=(mpfr_t& a);
  inline FP_NR<F>& operator=(mpfr_t &a) { set_mpfr(a, MPFR_RNDN); };

  inline bool operator==(const FP_NR<F> &a) const { return cmp(a) == 0; }
  inline bool operator==(double a) const { return cmp(a) == 0; }
  inline bool operator!=(const FP_NR<F> &a) const { return cmp(a) != 0; }
  inline bool operator!=(double a) const { return cmp(a) != 0; }

  inline bool operator<(const FP_NR<F> &a) const;
  inline bool operator<(double a) const;
  inline bool operator>(const FP_NR<F> &a) const;
  inline bool operator>(double a) const;
  inline bool operator<=(const FP_NR<F> &a) const;
  inline bool operator<=(double a) const;
  inline bool operator>=(const FP_NR<F> &a) const;
  inline bool operator>=(double a) const;

  inline FP_NR<F> &operator+=(const FP_NR<F> &a) { this->add(*this, a); return *this;};
  inline FP_NR<F> &operator-=(const FP_NR<F> &a) { this->sub(*this, a); return *this;};
  inline FP_NR<F> &operator*=(const FP_NR<F> &a) { this->mul(*this, a); return *this;};
  inline FP_NR<F> &operator/=(const FP_NR<F> &a) { this->div(*this, a); return *this;};

  inline FP_NR<F> &operator+=(const double a) { FP_NR<F> t=a; this->add(*this, t); return *this;};
  inline FP_NR<F> &operator-=(const double a) { FP_NR<F> t=a; this->sub(*this, t); return *this;};
  inline FP_NR<F> &operator*=(const double a) { this->mul_d(*this, a); return *this;};
  inline FP_NR<F> &operator/=(const double a) { this->div(*this, FP_NR<F>(a)); return *this;};


  /**
   * max between a and b
   */
  inline FP_NR &max_f(FP_NR<F> &b)
  {
    if ((*this) <= b)
      return b;
    else
      return (*this);
  }

  /**
   * Returns non-zero if the current value is zero, 0 otherwise.
   */
  inline bool is_zero() const;
  inline int zero_p() const { return is_zero(); }

  /**
   * Returns non-zero if the current value is NaN, 0 otherwise.
   */
  inline int is_nan() const;

  /**
   * Returns non-zero if !isnan(value) and !isinf(value), 0 otherwise.
   */
  inline int is_finite() const;

  /** arithmetic */

  /**
   * value := a + b.
   */
  inline void add(const FP_NR<F> &a, const FP_NR<F> &b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := a - b.
   */
  inline void sub(const FP_NR<F> &a, const FP_NR<F> &b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := a * b.
   */
  inline void mul(const FP_NR<F> &a, const FP_NR<F> &b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := a * b where b is double.
   */
  inline void mul_d(const FP_NR<F> &a, const double b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := a * b where b is mpfr_t.
   */
  inline void mul_mpfr(const FP_NR<mpfr_t> &a, const mpfr_t b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := a * 2^b.
   */
  inline void mul_2si(const FP_NR<F> &a, long b);

  /**
   * value := a / b.
   */
  inline void div(const FP_NR<F> &a, const FP_NR<F> &b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := a / b.
   */
  inline void div_d(const FP_NR<F> &a, const double b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := value + a * b.
   */
  inline void addmul(const FP_NR<F> &a, const FP_NR<F> &b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := value - a * b.
   */
  inline void submul(const FP_NR<F> &a, const FP_NR<F> &b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := a^b.
   */
  inline void pow_si(const FP_NR<F> &a, long b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := e^b.
   */
  inline void exponential(const FP_NR<F> &b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := natural logarithm of a.
   */
  inline void log(const FP_NR<F> &a, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := square root of b.
   */
  inline void sqrt(const FP_NR<F> &b, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := k-th root of b.
   */
  inline void root(const FP_NR<F> &b, unsigned int k, mp_rnd_t rnd = GMP_RNDN);

  /**
   * value := -b.
   */
  inline void neg(const FP_NR<F> &b);
  inline FP_NR<F> operator-() const
  {
    FP_NR<F> r;
    r.neg(*this);
    return r;
  }


  /**
   * value := absolute value of b.
   */
  inline void abs(const FP_NR<F> &b);

  /**
   * value := rounding of b to the nearest integer.
   */
  inline void rnd(const FP_NR<F> &b);

  /**
   * value <- (rounding of a * 2^expo_add) / 2^expo_add, but never overflows.
   * expo_add must be 0 if T=dpe_t or T=mpfr_t.
   */
  inline void rnd_we(const FP_NR<F> &b, long expo_add);

  /**
   * value := largest integer not greater than b.
   */
  inline void floor(const FP_NR<F> &b);

  /**
   * value := NaN.
   */
  inline void set_nan();

  /**
   * Efficiently swaps the values of two FP_NR.
   */
  inline void swap(FP_NR<F> &a);
  /**
   * value := sqrt(a^2 + b^2).
   * Needed for givens rotations
   */
  inline void hypot(const FP_NR<F> &a, const FP_NR<F> &b, mp_rnd_t rnd = GMP_RNDN);

};  // End class FP_NR

template <class F> inline FP_NR<F> operator+(const FP_NR<F> &a, const FP_NR<F> &b)
{
  FP_NR<F> r;
  r.add(a, b);
  return r;
}

template <class F> inline FP_NR<F>&& operator+(FP_NR<F> &&a, const FP_NR<F> &b)
{
  a.add(a, b);
  return std::move(a);
}


template <class F> inline FP_NR<F> operator+(const FP_NR<F> &a, double b)
{
  FP_NR<F> r;
  r = b;
  r.add(r, a);
  return r;
}

template <class F> inline FP_NR<F> operator+(double a, const FP_NR<F> &b)
{
  FP_NR<F> r;
  r = a;
  r.add(r, b);
  return r;
}

template <class F> inline FP_NR<F> operator+(const FP_NR<F> &a, long b)
{
  FP_NR<F> r;
  r = b;
  r.add(a, r);
  return r;
}

template <class F> inline FP_NR<F> operator+(long a, const FP_NR<F> &b)
{
  FP_NR<F> r;
  r = a;
  r.add(r, b);
  return r;
}


template <class F> inline FP_NR<F> operator-(const FP_NR<F> &a, const FP_NR<F> &b)
{
  FP_NR<F> r;
  r.sub(a, b);
  return r;
}


template <class F> inline FP_NR<F> operator-(const FP_NR<F> &a, double b)
{
  FP_NR<F> r;
  r = b;
  r.sub(a, r);
  return r;
}

template <class F> inline FP_NR<F> operator-(double a, const FP_NR<F> &b)
{
  FP_NR<F> r;
  r = a;
  r.sub(r, b);
  return r;
}

template <class F> inline FP_NR<F> operator-(const FP_NR<F> &a, long b)
{
  FP_NR<F> r;
  r = b;
  r.sub(a, r);
  return r;
}

template <class F> inline FP_NR<F> operator-(long a, const FP_NR<F> &b)
{
  FP_NR<F> r;
  r = a;
  r.sub(r, b);
  return r;
}

template <class F> inline FP_NR<F> operator*(const FP_NR<F> &a, const FP_NR<F> &b)
{
  FP_NR<F> r;
  r.mul(a, b);
  return r;
}

template <class F> inline FP_NR<F>&& operator*(FP_NR<F> &&a, const FP_NR<F> &b) {
  a.mul(a, b);
  return std::move(a);
}

template <class F> inline FP_NR<F> operator*(const FP_NR<F> &a, double b)
{
  FP_NR<F> r;
  r.mul_d(a, b);
  return r;
}

template <class F> inline FP_NR<F> operator*(const FP_NR<F> &a, long b)
{
  FP_NR<F> r;
  r = b;
  r.mul(r, a);
  return r;
}

template <class F> inline FP_NR<F> operator*(double a, const FP_NR<F> &b)
{
  FP_NR<F> r;
  r.mul_d(b, a);
  return r;
}


template <class F> inline FP_NR<F>&& operator*(FP_NR<F> &&a, double b) {
  a.mul_d(a, b);
  return std::move(a);
}

template <class F> inline FP_NR<F> operator/(const FP_NR<F> &a, const FP_NR<F> &b) {
  FP_NR<F> r;
  r.div(a, b);
  return r;
}

template <class F> inline FP_NR<F> operator/(const FP_NR<F> &a, const double b) {
  FP_NR<F> r;
  r.div(a, b);
  return r;
}

template <class F> inline FP_NR<F> operator/(const double a, const FP_NR<F> &b) {
  FP_NR<F> r;
  r.div(a, b);
  return r;
}

template <class F> inline FP_NR<F> sqrt(const FP_NR<F> &a) {
  FP_NR<F> r;
  r.sqrt(a);
  return r;
}

template <class F> inline FP_NR<F>&& sqrt(FP_NR<F> &&a) {
  a.sqrt(a);
  return std::move(a);
}

template <class F> inline FP_NR<F> abs(const FP_NR<F> &a) {
  FP_NR<F> r;
  r.abs(a);
  return r;
}

template <class F> inline FP_NR<F>&& abs(FP_NR<F> &&a) {
  a.abs(a);
  return std::move(a);
}


template <class F> inline FP_NR<F> floor(const FP_NR<F> &a) {
  FP_NR<F> r;
  r.floor(a);
  return r;
}

template <class F> inline FP_NR<F>&& floor(FP_NR<F> &&a) {
  a.floor(a);
  return std::move(a);
}


template <class F> inline FP_NR<F> log(const FP_NR<F> &a) {
  FP_NR<F> r;
  r.log(a);
  return r;
}

template <class F> inline FP_NR<F>&& log(FP_NR<F> &&a) {
  a.log(a);
  return std::move(a);
}

template <class F> inline FP_NR<F> exp(const FP_NR<F> &a) {
  FP_NR<F> r;
  r.exponential(a);
  return r;
}

template <class F> inline FP_NR<F>&& exp(FP_NR<F> &&a) {
  a.exponential(a);
  return std::move(a);
}

template <class F> inline FP_NR<F> pow_si(const FP_NR<F> &a, long b) {
  FP_NR<F> r;
  r.pow_si(a, b);
  return r;
}

template <class F> inline FP_NR<F> &&pow_si(FP_NR<F> &&a, long b) {
  a.pow_si(a, b);
  return std::move(a);
}

template <class F> inline FP_NR<F> root(const FP_NR<F> &a, unsigned int b) {
  FP_NR<F> r;
  r.root(a, b);
  return r;
}

template <class F> inline FP_NR<F> &&root(FP_NR<F> &&a, unsigned int b) {
  a.root(a, b);
  return std::move(a);
}


/**
 * FP_NR: some generic functions.
 */
template <class F> inline void FP_NR<F>::addmul(const FP_NR<F> &b, const FP_NR<F> &c, mp_rnd_t rnd)
{
  FP_NR<F> product;
  product.mul(b, c, rnd);
  add(*this, product, rnd);
}
template <class F> inline void FP_NR<F>::submul(const FP_NR<F> &b, const FP_NR<F> &c, mp_rnd_t rnd)
{
  FP_NR<F> product;
  product.mul(b, c, rnd);
  sub(*this, product, rnd);
}
template <class F> inline bool FP_NR<F>::operator<(const FP_NR<F> &a) const { return cmp(a) < 0; }
template <class F> inline bool FP_NR<F>::operator>(const FP_NR<F> &a) const { return cmp(a) > 0; }
template <class F> inline bool FP_NR<F>::operator<(double a) const { return cmp(a) < 0; }
template <class F> inline bool FP_NR<F>::operator>(double a) const { return cmp(a) > 0; }

/** overloading stream operators */

/**
 * Prints x on stream os.
 */
template <class T> ostream &operator<<(ostream &os, const FP_NR<T> &x) { return os << x.get_data(); }

#ifdef FPLLL_WITH_DPE
template <> ostream &operator<<(ostream &os, const FP_NR<dpe_t> &x);
#endif

template <> ostream &operator<<(ostream &os, const FP_NR<mpfr_t> &x);

FPLLL_END_NAMESPACE

#endif
