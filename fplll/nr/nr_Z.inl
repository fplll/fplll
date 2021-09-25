/*********************
 *   Integer Class
 *********************/

#ifndef FPLLL_NR_Z_H
#define FPLLL_NR_Z_H

FPLLL_BEGIN_NAMESPACE

/* declaration FP_NR */
template <class F> class FP_NR;

/**
 * Z_NR stores integers. This template provides a uniform interface
 *  for doing integer computations with several underlying types
 *  (long, double and mpz_t).
 */
template <class Z = mpz_t> class Z_NR
{

  Z data;

public:
  /**
   * Constructors.
   */
  inline Z_NR();
  inline Z_NR(const Z_NR<Z> &z);
  inline Z_NR(const Z &z);
  inline ~Z_NR();

  /** get data */

  /**
   * Returns the internal representation of the data.
   */
  inline Z &get_data() { return data; }

  /**
   * Returns the internal representation of the data.
   */
  inline const Z &get_data() const { return data; }

  /**
   * Converts this object to a double. If it does not fit in
   *  a double, the result is undefined.
   */
  inline double get_d() const;

  /**
   * Converts this object to a long double. If it does not fit in
   * a long double, the result is undefined.
   */
#ifdef FPLLL_WITH_LONG_DOUBLE
  inline long double get_ld() const;
#endif

  /**
   * Converts this object to a long. If it does not fit in a long,
   * the result is undefined.
   */
  inline long get_si() const;

  /**
   * Converts this object to a mpz.
   */
  inline void get_mpz(mpz_t r) const;

  /**
   * Returns the smallest non-negative expo such that |value| < 2^expo.
   */
  inline long exponent() const;

  /**
   * Computes f and expo such that 0.5 <= f < 1 and value ~= 2^expo * f.
   *   The rounding direction is undefined.
   *   This function is implemented only for native floating-point types
   *   (double and long double), otherwise the behaviour is undefined.
   */
  template <class F> inline void get_f_exp(F &f, long &expo);

  /** set data */

  /**
   * Sets the value to x. When FT=mpfr_t, x is rounded to the nearest
   * integer and if the fractional part of x is 0.5, the even integer
   * is chosen when. Otherwise, the rounding direction is undefined.
   */
  template <class F> inline void set_f(const FP_NR<F> &f);

  /**
   * Sets the value to s, signed integer in basis 10.
   */
  inline void set_str(const char *s);

  /** comparison */

  /**
   * 3-way comparison. Returns a positive number if *this > m, a
   * negative number if *this &lt; m or zero is *this == m.
   */
  inline int cmp(const Z_NR<Z> &m) const;

  /**
   * Sign. Returns a positive number, a negative number or zero
   * if the value of this object is respectively positive,
   * negative or null.
   */
  inline int sgn() const;

  /**
   * Operator
   */
  inline void operator=(const Z_NR<Z> &z);
  inline void operator=(const mpz_t &z);
  inline void operator=(long i);
  inline bool operator<(const Z_NR<Z> &a) const;
  inline bool operator<(long a) const;
  inline bool operator>(const Z_NR<Z> &a) const;
  inline bool operator>(long a) const;
  inline bool operator<=(const Z_NR<Z> &a) const;
  inline bool operator<=(long a) const;
  inline bool operator>=(const Z_NR<Z> &a) const;
  inline bool operator>=(long a) const;
  inline bool operator==(const Z_NR<Z> &a) const;
  inline bool operator==(long a) const;
  inline bool operator!=(const Z_NR<Z> &a) const;
  inline bool operator!=(long a) const;

  /**
   * max between a and b
   */
  inline Z_NR &max_z(Z_NR<Z> &b)
  {
    if ((*this) <= b)
      return b;
    else
      return (*this);
  }

  /**
   * Returns non-zero if the current value is zero, 0 otherwise.
   */
  inline bool is_zero() const { return *this == 0; }

  /** arithmetic */

  /**
   * value := a + b.
   */
  inline void add(const Z_NR<Z> &a, const Z_NR<Z> &b);
  inline void add_ui(const Z_NR<Z> &a, unsigned int b);
  /**
   * value := a - b.
   */
  inline void sub(const Z_NR<Z> &a, const Z_NR<Z> &b);
  inline void sub_ui(const Z_NR<Z> &a, unsigned int b);

  /**
   * value := a mod b
   */
  inline void mod(const Z_NR<Z> &a, const Z_NR<Z> &b);

  /**
   * value := -a.
   */
  inline void neg(const Z_NR<Z> &a);

  /**
   * value := a * b.
   */
  inline void mul(const Z_NR<Z> &a, const Z_NR<Z> &b);

  /**
   * value := a * b.
   */
  inline void mul_si(const Z_NR<Z> &a, long b);
  inline void mul_ui(const Z_NR<Z> &a, unsigned long b);

  /**
   * value := a * 2^b.
   * if Z=long and |b| >= size_in_bits(a), the result is undefined.
   */
  inline void mul_2si(const Z_NR<Z> &a, long b);

  /**
   * value := a / 2^b.
   * if Z=long and |b| >= size_in_bits(a), the result is undefined.
   */
  inline void div_2si(const Z_NR<Z> &a, long b);

  /**
   * value := value + a * b.
   */
  inline void addmul(const Z_NR<Z> &a, const Z_NR<Z> &b);
  inline void addmul_ui(const Z_NR<Z> &a, unsigned long b);
  inline void addmul_si(const Z_NR<Z> &a, long b);

  /**
   * value := value - a * b.
   */
  inline void submul(const Z_NR<Z> &a, const Z_NR<Z> &b);
  inline void submul_ui(const Z_NR<Z> &a, unsigned long b);

  /**
   * value := absolute value of a.
   */
  inline void abs(const Z_NR<Z> &a);

  /**
   * Efficiently swaps the values of two Z_NR.
   */
  inline void swap(Z_NR<Z> &a);

  /**
   * Generates random integer between 0 and 2^bits-1.
   */
  inline void randb(int bits);
  inline void randb_si(int bits);

  /**
   * Generates random integer between 0 and max - 1
   */
  inline void randm(const Z_NR<Z> &max);
  inline void randm_si(const Z_NR<Z> &max);

  /**
   * Generates smallest prime number above nbr
   */
  inline void nextprime(const Z_NR<Z> &nbr);
};

/**
 * LDConvHelper provides conversion functions between mpz_t and
 * long double which are not (yet) in GMP. It uses mpfr so it is slow.
 */
#ifdef FPLLL_WITH_LONG_DOUBLE

// Thread_local variables for doing numeric conversions.
// These are global to replace their previous usage, but member variables
// cannot be thread local and so they're declared at this scope.
// NOTE: these are extern because these can only have one definition.
// These are declared in util.cpp
extern thread_local mpfr_t temp_mpfr;
extern thread_local bool temp_mpfr_initialized;

class LDConvHelper
{
public:
  /** Converts op to a long double with rounding to nearest. */
  static long double mpz_get_ld(const mpz_t op)
  {
    init_temp();
    mpfr_set_z(temp_mpfr, op, GMP_RNDN);
    return mpfr_get_ld(temp_mpfr, GMP_RNDN);  // exact
  }

  static void free() { free_temp(); }

  /**
   * Returns d and sets exp such that 0.5 <= |d| < 1 and d * 2^exp is equal
   * to op rounded to the nearest long double.
   */
  static long double mpz_get_ld_2exp(long *exp, const mpz_t op)
  {
    init_temp();
    mpfr_set_z(temp_mpfr, op, GMP_RNDN);
    return mpfr_get_ld_2exp(exp, temp_mpfr, GMP_RNDN);  // exact
  }

  /** Sets the value of rop from op. */
  static void mpz_set_ld(mpz_t rop, long double op)
  {
    init_temp();
    mpfr_set_ld(temp_mpfr, op, GMP_RNDN);  // exact
    mpfr_get_z(rop, temp_mpfr, GMP_RNDN);
  }

private:
  static inline void init_temp()
  {
    if (!temp_mpfr_initialized)
    {
      mpfr_init2(temp_mpfr, numeric_limits<long double>::digits);
      temp_mpfr_initialized = true;
    }
  }

  static inline void free_temp()
  {
    if (temp_mpfr_initialized)
    {
      mpfr_clear(temp_mpfr);
      temp_mpfr_initialized = false;
    }
  }
};

#endif

/** overloading stream operators */

/**
 * Prints x on stream os.
 */
template <class T> ostream &operator<<(ostream &os, const Z_NR<T> &x);

/**
 * Reads x from stream is.
 */
template <class T> istream &operator>>(istream &is, Z_NR<T> &x);

FPLLL_END_NAMESPACE

#endif
