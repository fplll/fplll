/* Copyright (C) 2004, 2005, 2006, 2008 Patrick Pelissier, Paul Zimmermann,
   LORIA/INRIA Nancy - Grand-Est.

   This file is part of the DPE Library.

   The DPE Library is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

   The DPE Library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the DPE Library; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
   MA 02110-1301, USA. */

/* WARNING: Patched version */

#ifndef FPLLL_DPE_H
#define FPLLL_DPE_H
#define __DPE

#include <stdlib.h> /* For exit */
#include <stdio.h>  /* For fprintf and fscanf */
#include <math.h>   /* for round, floor, ceil */
#if defined (__sun) /* for round on Solaris 10 */ 
#include "tgmath.h"
#include "ieeefp.h"
#undef NAN
#define NAN 0.0/0.0
#endif 
#include <limits.h>

#define DPE_VERSION_MAJOR 1
#define DPE_VERSION_MINOR 5

#if defined(__GNUC__) && (__GNUC__ >= 3)
# define DPE_LIKELY(x) (__builtin_expect(!!(x),1))
# define DPE_UNLIKELY(x) (__builtin_expect((x),0))
# define DPE_UNUSED_ATTR  __attribute__((unused))
#else
# define DPE_LIKELY(x) (x)
# define DPE_UNLIKELY(x) (x)
# define DPE_UNUSED_ATTR
#endif

#if !defined(DPE_USE_DOUBLE) && !defined(DPE_USE_LONGDOUBLE)
# define DPE_USE_DOUBLE
#endif

#if (defined(__i386) || defined (__x86_64)) && !defined(DPE_LITTLEENDIAN32) && defined(DPE_USE_DOUBLE)
# define DPE_LITTLEENDIAN32
#endif

#if defined(DPE_USE_DOUBLE)
# define DPE_DOUBLE double /* mantissa type */
# define DPE_BITSIZE 53 /* bitsize of DPE_DOUBLE */
# define DPE_2_POW_BITSIZE 9007199254740992
/* DPE_LDEXP(DPE_DOUBLE m, DPEEXP e) return x = m * 2^e */
# define DPE_LDEXP ldexp
/* DPE_FREXP(DPE_DOUBLE x, DPEEXP *e) returns m, e such that x = m * 2^e with
   1/2 <= m < 1 */
# define DPE_FREXP frexp
/* DPE_ROUND(DPE_DOUBLE x) returns the nearest integer to x */
# define DPE_ROUND round
# define DPE_RINT rint
# define DPE_FLOOR floor
# define DPE_CEIL ceil
# define DPE_TRUNC trunc
#elif defined(DPE_USE_LONGDOUBLE)
# define DPE_DOUBLE long double
# define DPE_BITSIZE 64
# define DPE_2_POW_BITSIZE 18446744073709551616
# define DPE_LDEXP ldexpl
# define DPE_FREXP frexpl
# define DPE_ROUND roundl
# define DPE_RINT rintl
# define DPE_FLOOR floorl
# define DPE_CEIL ceill
# define DPE_TRUNC truncl
#else
# error "neither DPE_USE_DOUBLE nor DPE_USE_LONGDOUBLE is defined"
#endif

#if defined(DPE_USE_LONG)
# define DPE_EXP_T  long    /* exponent type */
# define DPE_EXPMIN LONG_MIN /* smallest possible exponent */
#elif defined(DPE_USE_LONGLONG)
# define DPE_EXP_T  long long
# define DPE_EXPMIN LONG_LONG_MIN
#else
# define DPE_EXP_T  int     /* exponent type */
# define DPE_EXPMIN INT_MIN /* smallest possible exponent */
#endif

typedef union
{
  double d;
  int i[2];
} dpe_double_words;

typedef struct
{
  DPE_DOUBLE d; /* significand */
  DPE_EXP_T exp; /* exponent */
} dpe_struct;

typedef dpe_struct dpe_t[1];

#define DPE_MANT(x) ((x)->d)
#define DPE_EXP(x)  ((x)->exp)
#define DPE_SIGN(x) ((DPE_MANT(x) < 0.0) ? -1 : (DPE_MANT(x) > 0.0))

#define DPE_INLINE static inline

/* initialize */
DPE_INLINE void
dpe_init (dpe_t x DPE_UNUSED_ATTR)
{
}

/* clear */
DPE_INLINE void
dpe_clear (dpe_t x DPE_UNUSED_ATTR)
{
}

/* set x to y */
DPE_INLINE void
dpe_set (dpe_t x, const dpe_t y)
{
  DPE_MANT(x) = DPE_MANT(y);
  DPE_EXP(x) = DPE_EXP(y);
}

/* set x to -y */
DPE_INLINE void
dpe_neg (dpe_t x, const dpe_t y)
{
  DPE_MANT(x) = -DPE_MANT(y);
  DPE_EXP(x) = DPE_EXP(y);
}

/* set x to |y| */
DPE_INLINE void
dpe_abs (dpe_t x, const dpe_t y)
{
  DPE_MANT(x) = (DPE_MANT(y) >= 0) ? DPE_MANT(y) : -DPE_MANT(y);
  DPE_EXP(x) = DPE_EXP(y);
}

/* set mantissa in [1/2, 1), except for 0 which has minimum exponent */
static void
dpe_normalize (dpe_t x)
{
  if (DPE_UNLIKELY (DPE_MANT(x) == 0.0 || std::isfinite (DPE_MANT(x)) == 0))
  {
    if (DPE_MANT(x) == 0.0)
      DPE_EXP(x) = DPE_EXPMIN;
    /* otherwise let the exponent of NaN, Inf unchanged */
  }
  else
  {
    DPE_EXP_T e;
#ifdef DPE_LITTLEENDIAN32 /* 32-bit little endian */
    dpe_double_words dw;
    dw.d = DPE_MANT(x);
    e = (dw.i[1] >> 20) & 0x7FF; /* unbiased exponent, 1022 for m=1/2 */
    DPE_EXP(x) += e - 1022;
    dw.i[1] = (dw.i[1] & 0x800FFFFF) | 0x3FE00000;
    DPE_MANT(x) = dw.d;
#else /* portable code */
    double m = DPE_MANT(x);
    DPE_MANT(x) = DPE_FREXP (m, &e);
    DPE_EXP(x) += e;
#endif
  }
}

#if defined(DPE_USE_DOUBLE)
static const double dpe_scale_tab[54] = {
1.0000000000000000000  , 0.5000000000000000000  , 0.2500000000000000000  , 
0.1250000000000000000  , 0.0625000000000000000  , 0.0312500000000000000  , 
0.0156250000000000000  , 0.0078125000000000000  , 0.0039062500000000000  , 
0.0019531250000000000  , 0.0009765625000000000  , 0.0004882812500000000  , 
0.0002441406250000000  , 0.0001220703125000000  , 0.0000610351562500000  , 
0.0000305175781250000  , 0.0000152587890625000  , 7.6293945312500000e-6  , 
3.8146972656250000e-6  , 1.9073486328125000e-6  , 9.5367431640625000e-7  , 
4.7683715820312500e-7  , 2.3841857910156250e-7  , 1.1920928955078125e-7  , 
5.9604644775390625e-8  , 2.9802322387695312e-8  , 1.4901161193847656e-8  , 
7.4505805969238281e-9  , 3.7252902984619141e-9  , 1.8626451492309570e-9  , 
9.3132257461547852e-10 , 4.6566128730773926e-10 , 2.3283064365386963e-10 , 
1.1641532182693481e-10 , 5.8207660913467407e-11 , 2.9103830456733704e-11 , 
1.4551915228366852e-11 , 7.2759576141834259e-12 , 3.6379788070917130e-12 , 
1.8189894035458565e-12 , 9.0949470177292824e-13 , 4.5474735088646412e-13 , 
2.2737367544323206e-13 , 1.1368683772161603e-13 , 5.6843418860808015e-14 , 
2.8421709430404007e-14 , 1.4210854715202004e-14 , 7.1054273576010019e-15 , 
3.5527136788005009e-15 , 1.7763568394002505e-15 , 8.8817841970012523e-16 , 
4.4408920985006262e-16 , 2.2204460492503131e-16 , 1.1102230246251565e-16 };
#endif

DPE_INLINE double
dpe_scale (double d, int s)
{
  /* -DPE_BITSIZE < s <= 0 and 1/2 <= d < 1 */
#if defined(DPE_USE_DOUBLE)
  return d * dpe_scale_tab [-s];
#else /* portable code */
  return DPE_LDEXP (d, s);
#endif
}

/* set x to y */
DPE_INLINE void
dpe_set_d (dpe_t x, double y)
{
  DPE_MANT(x) = (DPE_DOUBLE) y;
  DPE_EXP(x) = 0;
  dpe_normalize (x);
}

/* set x to y */
DPE_INLINE void
dpe_set_ld (dpe_t x, long double y)
{
  DPE_MANT(x) = (DPE_DOUBLE) y;
  DPE_EXP(x) = 0;
  dpe_normalize (x);
}

/* set x to y */
DPE_INLINE void
dpe_set_ui (dpe_t x, unsigned long y)
{
  DPE_MANT(x) = (DPE_DOUBLE) y;
  DPE_EXP(x) = 0;
  dpe_normalize (x);
}

/* set x to y */
DPE_INLINE void
dpe_set_si (dpe_t x, long y)
{
  DPE_MANT(x) = (DPE_DOUBLE) y;
  DPE_EXP(x) = 0;
  dpe_normalize (x);
}

DPE_INLINE long
dpe_get_si (const dpe_t x)
{
  DPE_DOUBLE d = DPE_LDEXP (DPE_MANT (x), DPE_EXP (x));
  return (long) d;
}

DPE_INLINE unsigned long
dpe_get_ui (const dpe_t x)
{
  DPE_DOUBLE d = DPE_LDEXP (DPE_MANT (x), DPE_EXP (x));
  return (d < 0.0) ? 0 : (unsigned long) d;
}

DPE_INLINE double
dpe_get_d (const dpe_t x)
{
  return DPE_LDEXP (DPE_MANT (x), DPE_EXP (x));
}

DPE_INLINE long double
dpe_get_ld (const dpe_t x)
{
  return DPE_LDEXP (DPE_MANT (x), DPE_EXP (x));
}


#ifdef __GMP_H__
/* set x to y */
DPE_INLINE void
dpe_set_z (dpe_t x, const mpz_t y)
{
  long e;
  DPE_MANT(x) = mpz_get_d_2exp (&e, y);
  DPE_EXP(x) = (DPE_EXP_T) e;
}

/* set x to y, rounded to nearest */
DPE_INLINE void
dpe_get_z (mpz_t x, const dpe_t y)
{
  DPE_EXP_T ey = DPE_EXP(y);
  if (ey >= DPE_BITSIZE) /* y is an integer */
  {
    DPE_DOUBLE d = DPE_MANT(y) * DPE_2_POW_BITSIZE; /* d is an integer */
    mpz_set_d (x, d); /* should be exact */
    mpz_mul_2exp (x, x, (unsigned long) ey - DPE_BITSIZE);
  }
  else /* DPE_EXP(y) < DPE_BITSIZE */
  {
    if (DPE_UNLIKELY (ey < 0)) /* |y| < 1/2 */
      mpz_set_ui (x, 0);
    else
    {
      DPE_DOUBLE d = DPE_LDEXP(DPE_MANT(y), ey);
      mpz_set_d (x, (double) DPE_ROUND(d));
    }
  }
}

/* return e and x such that y = x*2^e */
DPE_INLINE mp_exp_t
dpe_get_z_exp (mpz_t x, const dpe_t y)
{
  mpz_set_d (x, DPE_MANT (y) * DPE_2_POW_BITSIZE);
  return DPE_EXP(y) - DPE_BITSIZE;
}


#endif //#ifdef __GMP_H__


#ifdef __MPFR_H
/* set x to y where y is mpfr_t. */
DPE_INLINE void
dpe_get_mpfr (mpfr_t y, const dpe_t x, const mp_rnd_t rnd)
{
  mpfr_set_d (y, DPE_MANT(x), rnd);
  mpfr_mul_2si (y, y, DPE_EXP(x), rnd);
}

/* set y to x */
DPE_INLINE void
dpe_set_mpfr (dpe_t x, const mpfr_t y)
{
  long exp;
  DPE_MANT(x) = mpfr_get_d_2exp (&exp, y, GMP_RNDN);
  DPE_EXP(x) = exp;
  dpe_normalize (x);
}

#endif // __MPFR_H


/* x <- y + z, assuming y and z are normalized, returns x normalized */
DPE_INLINE void
dpe_add (dpe_t x, const dpe_t y, const dpe_t z)
{
  if (DPE_UNLIKELY (DPE_EXP(y) > DPE_EXP(z) + DPE_BITSIZE))
    /* |z| < 1/2*ulp(y), thus o(y+z) = y */
    dpe_set (x, y);
  else if (DPE_UNLIKELY (DPE_EXP(z) > DPE_EXP(y) + DPE_BITSIZE))
    dpe_set (x, z);
  else
  {
    DPE_EXP_T d = DPE_EXP(y) - DPE_EXP(z); /* |d| <= DPE_BITSIZE */

    if (d >= 0)
    {
      DPE_MANT(x) = DPE_MANT(y) + dpe_scale (DPE_MANT(z), -d);
      DPE_EXP(x) = DPE_EXP(y);
    }
    else
    {
      DPE_MANT(x) = DPE_MANT(z) + dpe_scale (DPE_MANT(y), d);
      DPE_EXP(x) = DPE_EXP(z);
    }
    dpe_normalize (x);
  }
}

/* x <- y - z, assuming y and z are normalized, returns x normalized */
DPE_INLINE void
dpe_sub (dpe_t x, const dpe_t y, const dpe_t z)
{
  if (DPE_UNLIKELY (DPE_EXP(y) > DPE_EXP(z) + DPE_BITSIZE))
    /* |z| < 1/2*ulp(y), thus o(y-z) = y */
    dpe_set (x, y);
  else if (DPE_UNLIKELY (DPE_EXP(z) > DPE_EXP(y) + DPE_BITSIZE))
    dpe_neg (x, z);
  else
  {
    DPE_EXP_T d = DPE_EXP(y) - DPE_EXP(z); /* |d| <= DPE_BITSIZE */

    if (d >= 0)
    {
      DPE_MANT(x) = DPE_MANT(y) - dpe_scale (DPE_MANT(z), -d);
      DPE_EXP(x) = DPE_EXP(y);
    }
    else
    {
      DPE_MANT(x) = dpe_scale (DPE_MANT(y), d) - DPE_MANT(z);
      DPE_EXP(x) = DPE_EXP(z);
    }
    dpe_normalize (x);
  }
}

/* x <- y * z, assuming y and z are normalized, returns x normalized */
DPE_INLINE void
dpe_mul (dpe_t x, const dpe_t y, const dpe_t z)
{
  DPE_MANT(x) = DPE_MANT(y) * DPE_MANT(z);
  DPE_EXP(x) = DPE_EXP(y) + DPE_EXP(z);
  dpe_normalize (x);
}

/* x <- y * z, assuming y and z are normalized, returns x normalized */
DPE_INLINE void
dpe_mul_d (dpe_t x, const dpe_t y, const double z)
{
  dpe_t zt;
  dpe_set_d (zt, z);
  DPE_MANT(x) = DPE_MANT(y) * DPE_MANT(zt);
  DPE_EXP(x) = DPE_EXP(y) + DPE_EXP(zt);
  dpe_normalize (x);
}

/* x <- x + y * z, assuming y and z are normalized, returns x normalized */
DPE_INLINE void
dpe_addmul (dpe_t x, const dpe_t y, const dpe_t z)
{
  dpe_t yz;
  dpe_mul (yz, y, z);
  dpe_add (x, x, yz);
}

/* x <- x - y * z, assuming y and z are normalized, returns x normalized */
DPE_INLINE void
dpe_submul (dpe_t x, const dpe_t y, const dpe_t z)
{
  dpe_t yz;
  dpe_mul (yz, y, z);
  dpe_sub (x, x, yz);
}

/* x <- sqrt(y), assuming y is normalized, returns x normalized */
DPE_INLINE void
dpe_sqrt (dpe_t x, const dpe_t y)
{
  DPE_EXP_T ey = DPE_EXP(y);
  if (ey % 2)
  {
    /* since 1/2 <= my < 1, 1/4 <= my/2 < 1 */
    DPE_MANT(x) = sqrt (0.5 * DPE_MANT(y));
    DPE_EXP(x) = (ey + 1) / 2;
  }
  else
  {
    DPE_MANT(x) = sqrt (DPE_MANT(y));
    DPE_EXP(x) = ey / 2;
  }
}

/* x <- y / z, assuming y and z are normalized, returns x normalized.
   Assumes z is not zero. */
DPE_INLINE void
dpe_div (dpe_t x, const dpe_t y, const dpe_t z)
{
  DPE_MANT(x) = DPE_MANT(y) / DPE_MANT(z);
  DPE_EXP(x) = DPE_EXP(y) - DPE_EXP(z);
  dpe_normalize (x);
}

/* x <- y * z, assuming y normalized, returns x normalized */
DPE_INLINE void
dpe_mul_ui (dpe_t x, const dpe_t y, unsigned long z)
{
  DPE_MANT(x) = DPE_MANT(y) * (DPE_DOUBLE) z;
  DPE_EXP(x) = DPE_EXP(y);
  dpe_normalize (x);
}

/* x <- y / z, assuming y normalized, z non-zero, returns x normalized */
DPE_INLINE void
dpe_div_ui (dpe_t x, const dpe_t y, unsigned long z)
{
  DPE_MANT(x) = DPE_MANT(y) / (DPE_DOUBLE) z;
  DPE_EXP(x) = DPE_EXP(y);
  dpe_normalize (x);
}

/* x <- y * 2^e */
DPE_INLINE void
dpe_mul_2si (dpe_t x, const dpe_t y, long e)
{
  DPE_MANT(x) = DPE_MANT(y);
  DPE_EXP(x) = DPE_EXP(y) + (DPE_EXP_T) e;
}

/* x <- y * 2^e */
DPE_INLINE void
dpe_mul_2exp (dpe_t x, const dpe_t y, unsigned long e)
{
  DPE_MANT(x) = DPE_MANT(y);
  DPE_EXP(x) = DPE_EXP(y) + (DPE_EXP_T) e;
}

/* x <- y / 2^e */
DPE_INLINE void
dpe_div_2exp (dpe_t x, const dpe_t y, unsigned long e)
{
  DPE_MANT(x) = DPE_MANT(y);
  DPE_EXP(x) = DPE_EXP(y) - (DPE_EXP_T) e;
}

/* return e and x such that y = x*2^e (equality is not guaranteed if the 'long'
   type has fewer bits than the significand in dpe_t) */
DPE_INLINE DPE_EXP_T
dpe_get_si_exp (long *x, const dpe_t y)
{
  if (sizeof(long) == 4) /* 32-bit word: long has 31 bits */
  {
    *x = (long) (DPE_MANT(y) * 2147483648.0);
    return DPE_EXP(y) - 31;
  }
  else if (sizeof(long) == 8) /* 64-bit word: long has 63 bits */
  {
    *x = (long) (DPE_MANT (y) * 9223372036854775808.0);
    return DPE_EXP(y) - 63;
  }
  else
  {
    fprintf (stderr, "Error, neither 32-bit nor 64-bit word\n");
    exit (1);
  }
}

static DPE_UNUSED_ATTR int dpe_str_prec = 16;
static int dpe_out_str (FILE *s, int base, const dpe_t x) DPE_UNUSED_ATTR;

static int
dpe_out_str (FILE *s, int base, const dpe_t x)
{
  DPE_DOUBLE d = DPE_MANT(x);
  DPE_EXP_T e2 = DPE_EXP(x);
  int e10 = 0;
  char sign = ' ';
  if (DPE_UNLIKELY (base != 10))
  {
    fprintf (stderr, "Error in dpe_out_str, only base 10 allowed\n");
    exit (1);
  }
  if (d == 0.0)
#ifdef DPE_USE_DOUBLE
    return fprintf (s, "%1.*f", dpe_str_prec, d);
#else
  return fprintf (s, "foo\n %1.*Lf", dpe_str_prec, d);
#endif
  if (d < 0)
  {
    d = -d;
    sign = '-';
  }
  if (e2 > 0)
  {
    while (e2 > 0)
    {
      e2 --;
      d *= 2.0;
      if (d >= 10.0)
      {
        d /= 10.0;
        e10 ++;
      }
    }
  }
  else /* e2 <= 0 */
  {
    while (e2 < 0)
    {
      e2 ++;
      d /= 2.0;
      if (d < 1.0)
      {
        d *= 10.0;
        e10 --;
      }
    }
  }
#ifdef DPE_USE_DOUBLE
  return fprintf (s, "%c%1.*f*10^%d", sign, dpe_str_prec, d, e10);
#else
  return fprintf (s, "%c%1.*Lf*10^%d", sign, dpe_str_prec, d, e10);
#endif
}

static size_t dpe_inp_str (dpe_t x, FILE *s, int base) DPE_UNUSED_ATTR;

static size_t
dpe_inp_str (dpe_t x, FILE *s, int base)
{
  size_t res;
  DPE_DOUBLE d;
  if (DPE_UNLIKELY (base != 10))
  {
    fprintf (stderr, "Error in dpe_out_str, only base 10 allowed\n");
    exit (1);
  }
#ifdef DPE_USE_DOUBLE
  res = fscanf (s, "%lf", &d);
#else
  res = fscanf (s, "%Lf", &d);
#endif
  dpe_set_d (x, d);
  return res;
}

DPE_INLINE void
dpe_dump (const dpe_t x)
{
  dpe_out_str (stdout, 10, x);
  putchar ('\n');
}

DPE_INLINE int
dpe_zero_p (const dpe_t x)
{
  return DPE_MANT (x) == 0;
}

/* return a positive value if x > y
   a negative value if x < y
   and 0 otherwise (x=y). */
DPE_INLINE int
dpe_cmp (const dpe_t x, const dpe_t y)
{
  int sx = DPE_SIGN(x);
  int d = sx - DPE_SIGN(y);

  if (d != 0)
    return d;
  else if (DPE_EXP(x) > DPE_EXP(y))
    return (sx > 0) ? 1 : -1;
  else if (DPE_EXP(y) > DPE_EXP(x))
    return (sx > 0) ? -1 : 1;
  else /* DPE_EXP(x) = DPE_EXP(y) */
    return (DPE_MANT(x) < DPE_MANT(y)) ? -1 : (DPE_MANT(x) > DPE_MANT(y));
}

DPE_INLINE int
dpe_cmp_d (const dpe_t x, double d)
{
  dpe_t y;
  dpe_set_d (y, d);
  return dpe_cmp (x, y);
}

DPE_INLINE int
dpe_cmp_ui (const dpe_t x, unsigned long d)
{
  dpe_t y;
  dpe_set_ui (y, d);
  return dpe_cmp (x, y);
}

DPE_INLINE int
dpe_cmp_si (const dpe_t x, long d)
{
  dpe_t y;
  dpe_set_si (y, d);
  return dpe_cmp (x, y);
}

/* set x to integer nearest to y */
DPE_INLINE void
dpe_round (dpe_t x, const dpe_t y)
{
  if (DPE_EXP(y) < 0) /* |y| < 1/2 */
    dpe_set_ui (x, 0);
  else if (DPE_EXP(y) >= DPE_BITSIZE) /* y is an integer */
    dpe_set (x, y);
  else
  {
    DPE_DOUBLE d;
    d = DPE_LDEXP(DPE_MANT(y), DPE_EXP(y));
    dpe_set_d (x, DPE_ROUND(d));
  }
}

/* set x to integer nearest to y */
DPE_INLINE void
dpe_rint (dpe_t x, const dpe_t y)
{
  if (DPE_EXP(y) < 0) /* |y| < 1/2 */
    dpe_set_ui (x, 0);
  else if (DPE_EXP(y) >= DPE_BITSIZE) /* y is an integer */
    dpe_set (x, y);
  else
  {
    DPE_DOUBLE d;
    d = DPE_LDEXP(DPE_MANT(y), DPE_EXP(y));
    dpe_set_d (x, DPE_RINT(d));
  }
}

/* set x to the fractional part of y, defined as y - trunc(y), thus the
   fractional part has absolute value in [0, 1), and same sign as y */
DPE_INLINE void
dpe_frac (dpe_t x, const dpe_t y)
{
  /* If |y| is smaller than 1, keep it */
  if (DPE_EXP(y) <= 0)
    dpe_set (x, y);
  else if (DPE_EXP(y) >= DPE_BITSIZE) /* y is an integer */
    dpe_set_ui (x, 0);
  else
  {
    DPE_DOUBLE d;
    d = DPE_LDEXP(DPE_MANT(y), DPE_EXP(y));
    dpe_set_d (x, d - DPE_TRUNC(d));
  }
}

/* set x to largest integer <= y */
DPE_INLINE void
dpe_floor (dpe_t x, const dpe_t y)
{
  if (DPE_EXP(y) <= 0) /* |y| < 1 */
  {
    if (DPE_SIGN(y) >= 0) /* 0 <= y < 1 */
      dpe_set_ui (x, 0);
    else /* -1 < y < 0 */
      dpe_set_si (x, -1);
  }
  else if (DPE_EXP(y) >= DPE_BITSIZE) /* y is an integer */
    dpe_set (x, y);
  else
  {
    DPE_DOUBLE d;
    d = DPE_LDEXP(DPE_MANT(y), DPE_EXP(y));
    dpe_set_d (x, DPE_FLOOR(d));
  }
}

/* set x to smallest integer >= y */
DPE_INLINE void
dpe_ceil (dpe_t x, const dpe_t y)
{
  if (DPE_EXP(y) <= 0) /* |y| < 1 */
  {
    if (DPE_SIGN(y) > 0) /* 0 < y < 1 */
      dpe_set_ui (x, 1);
    else /* -1 < y <= 0 */
      dpe_set_si (x, 0);
  }
  else if (DPE_EXP(y) >= DPE_BITSIZE) /* y is an integer */
    dpe_set (x, y);
  else
  {
    DPE_DOUBLE d;
    d = DPE_LDEXP(DPE_MANT(y), DPE_EXP(y));
    dpe_set_d (x, DPE_CEIL(d));
  }
}

DPE_INLINE void
dpe_swap (dpe_t x, dpe_t y)
{
  DPE_EXP_T i = DPE_EXP (x);
  DPE_DOUBLE d = DPE_MANT (x);
  DPE_EXP (x) = DPE_EXP (y);
  DPE_MANT (x) = DPE_MANT (y);
  DPE_EXP (y) = i;
  DPE_MANT (y) = d;
}

/* Ugly hacks: No correct rounding guaranteed */

DPE_INLINE void
dpe_ugly_log (dpe_t x, const dpe_t y)
{
  dpe_set_d (x, ((double) DPE_EXP(y)) *  M_LN2 + log (DPE_MANT(y)));
}

DPE_INLINE void
dpe_ugly_exp (dpe_t x, const dpe_t y)
{
  //printf ("## exp is %ld\n", DPE_EXP(y));
  dpe_set_d (x, exp(((double) DPE_MANT(y)) * pow(2, ((double)DPE_EXP(y))))); 
}


/* More hacks */
/* x = y^k */
DPE_INLINE void
dpe_pow_si (dpe_t x, const dpe_t y, const unsigned int k)
{
  DPE_MANT (x) = pow(DPE_MANT (y), k);
  DPE_EXP (x) = DPE_EXP (y) * k;
  dpe_normalize (x);
}

#ifdef __MPFR_H
/* x = e^y */
DPE_INLINE void
dpe_exponential (dpe_t x, const dpe_t y)
{
  /* floor(log(DBL_MAX)) = 709 */
  if (dpe_cmp_ui(y, 709) <= 0 && dpe_cmp_si(y, -709) >= 0)
    dpe_ugly_exp (x, y);
  else {
    mpfr_t t, s;
    mpfr_init (t);
    mpfr_init (s);
    mpfr_set_d (t, 2.0, GMP_RNDN);
    mpfr_pow_si (t, t, DPE_EXP(y), GMP_RNDN);
    mpfr_set_d (s, DPE_MANT(y), GMP_RNDN);
    mpfr_mul (t, t, s, GMP_RNDN);
    mpfr_exp (t, t, GMP_RNDN);
    dpe_set_mpfr (x, t);
    mpfr_clear (t);
    mpfr_clear (s);
    mpfr_free_cache();
  }
}
#endif

DPE_INLINE void
dpe_log (dpe_t x, const dpe_t y)
{
  dpe_set_d (x, (double) DPE_EXP(y) *  log(2.0) + log (DPE_MANT(y)));
}



#endif /* FPLLL_DPE_H */
