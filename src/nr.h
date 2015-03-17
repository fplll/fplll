/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2011 Xavier Pujol.

   This file is part of fplll. fplll is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   fplll is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with fplll. If not, see <http://www.gnu.org/licenses/>. */

#ifndef FPLLL_NR_H
#define FPLLL_NR_H

#include "defs.h"

FPLLL_BEGIN_NAMESPACE

template<class F> class FP_NR;

#ifdef FPLLL_WITH_LONG_DOUBLE

/**
 * LDConvHelper provides conversion functions between mpz_t and long double
 * which are not (yet) in GMP. It uses mpfr so it is slow.
 */
class LDConvHelper {
public:
  /** Converts op to a long double with rounding to nearest. */
  static long double mpz_get_ld(const mpz_t op) {
    initTemp();
    mpfr_set_z(temp, op, GMP_RNDN);
    return mpfr_get_ld(temp, GMP_RNDN); // exact
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
    if (!tempInitialized) {
      mpfr_init2(temp, numeric_limits<long double>::digits);
      tempInitialized = true;
    }
  }

  // These static members are initialized in util.cpp
  static mpfr_t temp;
  static bool tempInitialized;
};

#endif

/* Integers
   ======== */

/** Z_NR stores integers. This template provides a uniform interface for doing
    integer computations with several underlying types (long, double and
    mpz_t). */
template<class Z>
class Z_NR
{
  Z data;
public:
  /** Default constructor. The initial value is undefined. */
  inline Z_NR<Z>();
  /** Copy constructor. */
  inline Z_NR<Z>(const Z_NR<Z>& z);
  /** Destructor. */
  inline ~Z_NR<Z>();

  /** Converts this object to a double. If it does not fit in a double, the
      result is undefined. */
  inline double get_d() const;
#ifdef FPLLL_WITH_LONG_DOUBLE
  /** Converts this object to a long double. If it does not fit in a long
      double, the result is undefined. */
  inline long double get_ld() const;
#endif
  /** Converts this object to a long. If it does not fit in a long, the result
      is undefined. */
  inline long get_si() const;
  /** Computes f and expo such that 0.5 <= f < 1 and value ~= 2^expo * f.
      The rounding direction is undefined.
      This function is implemented only for native floating-point types
      (double and long double), otherwise the behaviour is undefined. */
  template<class F> inline void get_f_exp(F& f, long& expo);
  /** Sets the value to x. When FT=mpfr_t, x is rounded to the nearest integer
      and if the fractional part of x is 0.5, the even integer is chosen when.
      Otherwise, the rounding direction is undefined. */
  template<class F> inline void set_f(const FP_NR<F>& f);
  /** Sets the value to s, signed integer in basis 10. */
  inline void set_str(const char* s);
  /** Sets the value to z. */
  inline void operator=(const Z_NR<Z>& z);
  /** value := i. */
  inline void operator=(long i);

  /** 3-way comparison. Returns a positive number if *this > m, a negative
      number if *this &lt; m or zero is *this == m. */
  inline int cmp(const Z_NR<Z>& m) const;
  /** Sign. Returns a positive number, a negative number or zero if the value of
      this object is respectively positive, negative or null. */
  inline int sgn() const;
  /** Comparison operator. */
  inline bool operator<(const Z_NR<Z>& a) const;
  /** Comparison operator. */
  inline bool operator<(long a) const;
  /** Comparison operator. */
  inline bool operator>(const Z_NR<Z>& a) const;
  /** Comparison operator. */
  inline bool operator>(long a) const;
  /** Comparison operator. */
  inline bool operator<=(const Z_NR<Z>& a) const;
  /** Comparison operator. */
  inline bool operator<=(long a) const;
  /** Comparison operator. */
  inline bool operator>=(const Z_NR<Z>& a) const;
  /** Comparison operator. */
  inline bool operator>=(long a) const;
  /** Comparison operator. */
  inline bool operator==(const Z_NR<Z>& a) const;
  /** Comparison operator. */
  inline bool operator==(long a) const;
  /** Comparison operator. */
  inline bool operator!=(const Z_NR<Z>& a) const;
  /** Comparison operator. */
  inline bool operator!=(long a) const;

  /** value := a + b. */
  inline void add(const Z_NR<Z>& a, const Z_NR<Z>& b);
  inline void add_ui(const Z_NR<Z>& a, unsigned int b);
  /** value := a - b. */
  inline void sub(const Z_NR<Z>& a, const Z_NR<Z>& b);
  /** value := -a. */
  inline void neg(const Z_NR<Z>& a);
  /** value := a * b. */
  inline void mul(const Z_NR<Z>& a, const Z_NR<Z>& b);
  /** value := a * b. */
  inline void mul_si(const Z_NR<Z>& a, long b);
  inline void mul_ui(const Z_NR<Z>& a, unsigned long b);
  /** value := a * 2^b.
      If Z=long and |b| >= size_in_bits(a), the result is undefined. */
  inline void mul_2si(const Z_NR<Z>& a, long b);
  /** value := a / 2^b.
      If Z=long and |b| >= size_in_bits(a), the result is undefined. */
  inline void div_2si(const Z_NR<Z>& a, long b);

  /** value := value + a * b. */
  inline void addmul(const Z_NR<Z>& a, const Z_NR<Z>& b);
  inline void addmul_ui(const Z_NR<Z>& a, unsigned long b);
  inline void addmul_si(const Z_NR<Z>& a, long b);
  /** value := value - a * b. */
  inline void submul(const Z_NR<Z>& a, const Z_NR<Z>& b);
  inline void submul_ui(const Z_NR<Z>& a, unsigned long b);

  /** value := absolute value of a. */
  inline void abs(const Z_NR<Z>& a);
  /** Returns the smallest non-negative expo such that |value| < 2^expo. */
  inline long exponent() const;
  
  /** Generates random integer between 0 and 2^n-1. */
  inline void randb(int bits);
  /** Generates random integer between 0 and n - 1  */
  inline void randm(const Z_NR<Z>& max);

  inline bool is_zero() const {return *this == 0;}

  /** Efficiently swaps the values of two Z_NR. */
  inline void swap(Z_NR<Z>& a);

  /** Returns the internal representation of the data. */
  inline Z& getData()             {return data;}
  /** Returns the internal representation of the data. */
  inline const Z& getData() const {return data;}

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

/* Floats
   ====== */

/** FP_NR stores floating-point numbers. This template provides a uniform
    interface for doing floating-point computations with several underlying
    types (double, dpe_t and mpfr_t). For all functions, the rounding mode rnd
    is ignored unless F=mpfr_t. */
template<class F>
class FP_NR
{
  F data;
public:
  /** Default constructor. The initial value is undefined. */
  inline FP_NR<F>();
  /** Copy constructor. */
  inline FP_NR<F>(const FP_NR<F>& f);
  /** Destructor. */
  inline ~FP_NR<F>();

  /** Converts this object to a double. If it does not fit in a double, the
      result is undefined. */
  inline double get_d(mp_rnd_t rnd = GMP_RNDN) const;

  /** Converts this object to a long. The rounding direction is undefined.
      If it does not fit in a long, the result is undefined. */
  inline long get_si() const;

  /** Returns x and defines expo such that trunc(value) ~= x * 2^expo.
      The '~=' is an equality if trunc(value) <= LONG_MAX.
      expo is the minimum non-negative value such that |x| <= LONG_MAX. */
  inline long get_si_exp(long& expo) const;

  /** Returns x and defines expo such that trunc(value * 2^expoAdd) ~= x * 2^expo
      The '~=' is an equality if trunc(value * 2^expoAdd) <= LONG_MAX.
      expo is the minimum non-negative value such that x <= LONG_MAX.
      expoAdd must be 0 if T=dpe_t or T=mpfr_t. */
  inline long get_si_exp_we(long& expo, long expoAdd) const;

  /** Computes a and expo such that trunc(value) ~= a * 2^expo.
      The '~=' is an equality if Z=mpz_t. expo is always non-negative.
      Note that expo might be very close to LONG_MAX when value = 0. */
  template<class Z> inline void get_z_exp(Z_NR<Z>& a, long& expo) const;

  /** Computes a and expo such that trunc(value) * 2^expoAdd ~= a * 2^expo.
      The '~=' is an equality if Z=mpz_t. expo is always non-negative.
      expoAdd must be 0 if T=dpe_t or T=mpfr_t. */
  template<class Z> inline void get_z_exp_we(Z_NR<Z>& a, long& expo, long expoAdd) const;

  /** Sets the value to z. */
  template<class Z> inline void set_z(const Z_NR<Z>& z, mp_rnd_t rnd = GMP_RNDN);
  /** Sets the value to a. */
  inline void operator=(const FP_NR<F>& a);
  /** Sets the value to a. */
  inline void operator=(double a);

  /** 3-way comparison. Returns a positive number if *this > b, a negative
      number if *this &lt; b or zero is *this == b. */
  inline int cmp(const FP_NR<F>& b) const;
  /** 3-way comparison. Returns a positive number if *this > d, a negative
      number if *this &lt; d or zero is *this == d. */
  inline int cmp(double d) const;
  /** Sign. Returns a positive number, a negative number or zero if the value of
      this object is respectively positive, negative or null. */
  inline int sgn() const;
  /** Comparison operator. */
  inline bool operator<(const FP_NR<F>& a) const;
  /** Comparison operator. */
  inline bool operator<(double a) const;
  /** Comparison operator. */
  inline bool operator>(const FP_NR<F>& a) const;
  /** Comparison operator. */
  inline bool operator>(double a) const;
  /** Comparison operator. */
  inline bool operator<=(const FP_NR<F>& a) const;
  /** Comparison operator. */
  inline bool operator<=(double a) const;
  /** Comparison operator. */
  inline bool operator>=(const FP_NR<F>& a) const;
  /** Comparison operator. */
  inline bool operator>=(double a) const;

  /** value := a + b. */
  inline void add(const FP_NR<F>& b, const FP_NR<F>& c, mp_rnd_t rnd = GMP_RNDN);
  /** value := a - b. */
  inline void sub(const FP_NR<F>& b, const FP_NR<F>& c, mp_rnd_t rnd = GMP_RNDN);
  /** value := -b. */
  inline void neg(const FP_NR<F>& b);
  /** value := a * b. */
  inline void mul(const FP_NR<F>& b, const FP_NR<F>& c, mp_rnd_t rnd = GMP_RNDN);
  /** value := a * 2^b. */
  inline void mul_2si(const FP_NR<F>& b, long c);
  /** value := a / b. */
  inline void div(const FP_NR<F>& b, const FP_NR<F>& c, mp_rnd_t rnd = GMP_RNDN);
  /** value := value - a * b. */
  inline void submul(const FP_NR<F>& b, const FP_NR<F>& c, mp_rnd_t rnd = GMP_RNDN);
  /** value := value + a * b. */
  inline void addmul(const FP_NR<F>& b, const FP_NR<F>& c, mp_rnd_t rnd = GMP_RNDN);

  /** value := absolute value of b. */
  inline void abs(const FP_NR<F>& b);
  /** value := rounding of b to the nearest integer. */
  inline void rnd(const FP_NR<F>& b);
  /** value <- (rounding of a * 2^expoAdd) / 2^expoAdd, but never overflows.
      expoAdd must be 0 if T=dpe_t or T=mpfr_t. */
  inline void rnd_we(const FP_NR<F>& b, long expoAdd);
  /** value := largest integer not greater than b. */
  inline void floor(const FP_NR<F>& b);
  /** Returns expo such value 2^(expo-1) <= value < 2^expo
      (expo = floor(log2(value)) - 1).
      The return value is undefined if *this == 0. */
  inline long exponent() const;
  /** Returns non-zero if the current value is zero, 0 otherwise. */
  inline bool is_zero() const;
  inline int zero_p() const {return is_zero();}
  /** value := NaN. */
  inline void set_nan();
  /** Returns non-zero if the current value is NaN, 0 otherwise. */
  inline int is_nan() const;
  /** Returns non-zero if !isnan(value) and !isinf(value), 0 otherwise */
  inline int is_finite() const;
  /** value := square root of b. */
  inline void sqrt(const FP_NR<F>& b, mp_rnd_t rnd = GMP_RNDN);
  /** value := a * 2^b. */
  inline void pow_si(const FP_NR<F>& a, long b, mp_rnd_t rnd = GMP_RNDN);
  /** value := exponential of a. */
  inline void exponential(const FP_NR<F>& a, mp_rnd_t rnd = GMP_RNDN);
  /** value := natural logarithm of a. */
  inline void log(const FP_NR<F>& a, mp_rnd_t rnd = GMP_RNDN);
  /** Efficiently swaps the values of two FP_NR. */
  inline void swap(FP_NR<F>& a);

  /** Sets the precision of new FP_NR&lt;F&gt; objects. Returns the previous
      value. This has no effect is F != mpfr_t. */
  static inline unsigned int setprec(unsigned int prec);
  /** Returns the current precision for new FP_NR&lt;F&gt; objects. */
  static inline unsigned int getprec();

  /** Returns the internal representation of the data. */
  inline F& getData()             {return data;}
  /** Returns the internal representation of the data. */
  inline const F& getData() const {return data;}

#ifdef FPLLL_V3_COMPAT
  // Old interface (do not use)
  inline void print() const;
  inline void printerr() const;
  inline double get() const          {return get_d();}
  inline void set(const FP_NR<F>& s) {*this = s;}
  inline void set(double s)          {*this = s;}
  inline void set(unsigned int s);
  inline void mul_2ui(const FP_NR<F>& b, unsigned int c) {
    mul_2si(b, static_cast<long>(c));
  }
  inline void div_2ui(const FP_NR<F>& b, unsigned int c) {
    mul_2si(b, -static_cast<long>(c));
  }
  inline int exp() const          {return static_cast<int>(exponent());};
  inline F& GetData()             {return data;}
  inline const F& GetData() const {return data;}
#endif
};

/* Random generator for mpz_t
   ========================== */

class RandGen {
public:
  static void init() {
    initialized = true;
    gmp_randinit_default(gmpState);
  }
  static void initWithSeed(unsigned long seed) {
    init();
    gmp_randseed_ui(gmpState, seed);
  }
  static void initWithTime() {
    init();
    gmp_randseed_ui(gmpState, time(NULL));
  }
  static bool getInitialized() {
    return initialized;
  }
  static gmp_randstate_t& getGMPState() {
    if (!initialized) init();
    return gmpState;
  }
private:
  static bool initialized;
  static gmp_randstate_t gmpState;
};

template<class T> inline const char* numTypeStr()       {return "";}
#ifdef FPLLL_WITH_ZLONG
template<> inline const char* numTypeStr<long>()        {return "long";}
#endif
template<> inline const char* numTypeStr<double>()      {return "double";}
template<> inline const char* numTypeStr<mpz_t>()       {return "mpz_t";}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<> inline const char* numTypeStr<long double>() {return "long double";}
#endif
#ifdef FPLLL_WITH_DPE
template<> inline const char* numTypeStr<dpe_t>()       {return "dpe_t";}
#endif
template<> inline const char* numTypeStr<mpfr_t>()      {return "mpfr_t";}

typedef mpfr_t FloatT;

/* Default floating-point type */
typedef FP_NR<FloatT> Float;

typedef mpz_t IntegerT;

/* Default integer type */
typedef Z_NR<IntegerT> Integer;

/* Floating-point type inside the SVP/CVP solver */
typedef double enumf;

FPLLL_END_NAMESPACE

#include "nr.cpp"

#endif
