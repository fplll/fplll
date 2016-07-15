/*********************
 *   Integer Class
 *********************/

#ifndef FPLLL_NR_Z_H
#define FPLLL_NR_Z_H

FPLLL_BEGIN_NAMESPACE

/* declaration FP_NR */
template<class F> class FP_NR;


/**
 * Z_NR stores integers. This template provides a uniform interface 
 *  for doing integer computations with several underlying types
 *  (long, double and mpz_t). 
 */
template<class Z>
class Z_NR
{

  Z data;

public:

  /**
   * Constructors.
   */
  inline Z_NR<Z>();
  inline Z_NR<Z>(const Z_NR<Z>& z);
  inline ~Z_NR<Z>();

  /** get data */

  /**
   * Returns the internal representation of the data.
   */
  inline Z& get_data()             {return data;}

  /**
   * Returns the internal representation of the data.
   */
  inline const Z& get_data() const {return data;}

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
  template<class F> inline void get_f_exp(F& f, long& expo);

  /** set data */

  /**
   * Sets the value to x. When FT=mpfr_t, x is rounded to the nearest 
   * integer and if the fractional part of x is 0.5, the even integer 
   * is chosen when. Otherwise, the rounding direction is undefined.
   */
  template<class F> inline void set_f(const FP_NR<F>& f);

  /**
   * Sets the value to s, signed integer in basis 10. 
   */
  inline void set_str(const char* s);

  /** comparison */

  /**
   * 3-way comparison. Returns a positive number if *this > m, a 
   * negative number if *this &lt; m or zero is *this == m. 
   */
  inline int cmp(const Z_NR<Z>& m) const;

  /**
   * Sign. Returns a positive number, a negative number or zero 
   * if the value of this object is respectively positive, 
   * negative or null. 
   */
  inline int sgn() const;

  /**
   * Operator
   */
  inline void operator=(const Z_NR<Z>& z);
  inline void operator=(const mpz_t& z);
  inline void operator=(long i);
  inline bool operator<(const Z_NR<Z>& a) const;
  inline bool operator<(long a) const;
  inline bool operator>(const Z_NR<Z>& a) const;
  inline bool operator>(long a) const;
  inline bool operator<=(const Z_NR<Z>& a) const;
  inline bool operator<=(long a) const;
  inline bool operator>=(const Z_NR<Z>& a) const;
  inline bool operator>=(long a) const;
  inline bool operator==(const Z_NR<Z>& a) const;
  inline bool operator==(long a) const;
  inline bool operator!=(const Z_NR<Z>& a) const;
  inline bool operator!=(long a) const;

  /**
   * max between a and b
   */
  inline Z_NR& max_z(Z_NR<Z>& b) {
    if ((*this)<=b)
      return b;
    else
      return (*this);
  }

  /**
   * Returns non-zero if the current value is zero, 0 otherwise.
   */
  inline bool is_zero() const {return *this == 0;}

  /** arithmetic */

  /**
   * value := a + b. 
   */
  inline void add(const Z_NR<Z>& a, const Z_NR<Z>& b);
  inline void add_ui(const Z_NR<Z>& a, unsigned int b);
  /**
   * value := a - b. 
   */
  inline void sub(const Z_NR<Z>& a, const Z_NR<Z>& b);
  inline void sub_ui(const Z_NR<Z>& a, unsigned int b);

  /**
   * value := -a. 
   */
  inline void neg(const Z_NR<Z>& a);

  /**
   * value := a * b. 
   */
  inline void mul(const Z_NR<Z>& a, const Z_NR<Z>& b);

  /**
   * value := a * b. 
   */
  inline void mul_si(const Z_NR<Z>& a, long b);
  inline void mul_ui(const Z_NR<Z>& a, unsigned long b);

  /**
   * value := a * 2^b.
   * if Z=long and |b| >= size_in_bits(a), the result is undefined.
   */
  inline void mul_2si(const Z_NR<Z>& a, long b);

  /**
   * value := a / 2^b.
   * if Z=long and |b| >= size_in_bits(a), the result is undefined.
   */
  inline void div_2si(const Z_NR<Z>& a, long b);

  /**
   * value := value + a * b. 
   */
  inline void addmul(const Z_NR<Z>& a, const Z_NR<Z>& b);
  inline void addmul_ui(const Z_NR<Z>& a, unsigned long b);
  inline void addmul_si(const Z_NR<Z>& a, long b);

  /**
   * value := value - a * b. 
   */
  inline void submul(const Z_NR<Z>& a, const Z_NR<Z>& b);
  inline void submul_ui(const Z_NR<Z>& a, unsigned long b);

  /**
   * value := absolute value of a. 
   */
  inline void abs(const Z_NR<Z>& a);

  /**
   * Efficiently swaps the values of two Z_NR.
   */
  inline void swap(Z_NR<Z>& a);

  /**
   * Generates random integer between 0 and 2^bits-1. 
   */
  inline void randb(int bits);
  inline void randb_si(int bits);

  /**
   * Generates random integer between 0 and max - 1  
   */
  inline void randm(const Z_NR<Z>& max);
  inline void randm_si(const Z_NR<Z>& max);

  /**	 
   * Generates smallest prime number above nbr 
  */
  inline void nextprime(const Z_NR<Z>& nbr);

#ifdef FPLLL_V3_COMPAT
  // Old interface (do not use)
  inline void print() const;
  inline void printerr() const;
  inline void read();
  inline double get_d_2exp(long* expo) const;
  inline void set_si(long s) {*this = s;}
  inline void set(const Z_NR<Z>& s)     {*this = s;}
  inline void set(/*const*/ Z& s);
  inline void set(unsigned long s);
  inline void mul_2exp(const Z_NR<Z>& a, long expo) {mul_2si(a, expo);}
  inline void div_2exp(const Z_NR<Z>& a, long expo) {div_2si(a, expo);}
  inline Z& GetData()             {return data;}
  inline const Z& GetData() const {return data;}
#endif

};


/**
 * LDConvHelper provides conversion functions between mpz_t and
 * long double which are not (yet) in GMP. It uses mpfr so it is slow.
 */
#ifdef FPLLL_WITH_LONG_DOUBLE

class LDConvHelper {
public:
  /** Converts op to a long double with rounding to nearest. */
  static long double mpz_get_ld(const mpz_t op) {
    initTemp();
    mpfr_set_z(temp, op, GMP_RNDN);
    return mpfr_get_ld(temp, GMP_RNDN); // exact
  }

  static void free() {
    freeTemp();
  }

  /**
   * Returns d and sets exp such that 0.5 <= |d| < 1 and d * 2^exp is equal
   * to op rounded to the nearest long double.
   */
  static long double mpz_get_ld_2exp(long* exp, const mpz_t op) {
    initTemp();
    mpfr_set_z(temp, op, GMP_RNDN);
    return mpfr_get_ld_2exp(exp, temp, GMP_RNDN); // exact
  }

  /** Sets the value of rop from op. */
  static void mpz_set_ld(mpz_t rop, long double op) {
    initTemp();
    mpfr_set_ld(temp, op, GMP_RNDN); // exact
    mpfr_get_z(rop, temp, GMP_RNDN);
  }

private:
  static inline void initTemp() {
    if (!temp_initialized) {
      mpfr_init2(temp, numeric_limits<long double>::digits);
      temp_initialized = true;
    }
  }

  static inline void freeTemp() {
    if (temp_initialized) {
      mpfr_clear(temp);
      temp_initialized = false;
    }
  }

  // These static members are initialized in util.cpp
  static mpfr_t temp;
  static bool temp_initialized;
};

#endif


/** overloading stream operators */

/**
 * Prints x on stream os. 
 */
template<class T> 
ostream& operator<<(ostream& os, const Z_NR<T>& x);

/**
 * Reads x from stream is.
 */
template<class T>
istream& operator>>(istream& is, Z_NR<T>& x);


FPLLL_END_NAMESPACE

#endif
