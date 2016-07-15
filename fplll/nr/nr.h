#ifndef FPLLL_NR_H
#define FPLLL_NR_H

/**
 * loading definitions
 */
#include <iostream>
#include "../defs.h"

#include "nr_rand.inl"

#include "nr_Z.inl"
#include "nr_Z_l.inl"
#include "nr_Z_d.inl"
#include "nr_Z_mpz.inl"

#include "nr_FP.inl"
#include "nr_FP_d.inl"
#include "nr_FP_ld.inl"

#ifdef FPLLL_WITH_DPE
#include "nr_FP_dpe.inl"
#endif

#ifdef FPLLL_WITH_QD
#include "nr_FP_dd.inl"
#include "nr_FP_qd.inl"
#endif

#include "nr_FP_mpfr.inl"

#include "nr_Z_misc.inl"
#include "nr_FP_misc.inl"

FPLLL_BEGIN_NAMESPACE


/**
 * Some defaults.
 */

/* Default floating-point type */
typedef mpfr_t FloatT;
typedef FP_NR<FloatT> Float;

/* Default integer type */
typedef mpz_t IntegerT;
typedef Z_NR<IntegerT> Integer;

/* Floating-point type inside the SVP/CVP solver */
typedef double enumf;
typedef double enumxt;
//typedef int enumxt;

/**
 * return type
 */
template<class T> inline const char* num_type_str()       {return "";}
#ifdef FPLLL_WITH_ZLONG
template<> inline const char* num_type_str<long>()        {return "long";}
#endif
template<> inline const char* num_type_str<double>()      {return "double";}
template<> inline const char* num_type_str<mpz_t>()       {return "mpz_t";}
#ifdef FPLLL_WITH_LONG_DOUBLE
template<> inline const char* num_type_str<long double>() {return "long double";}
#endif
#ifdef FPLLL_WITH_DPE
template<> inline const char* num_type_str<dpe_t>()       {return "dpe_t";}
#endif
#ifdef FPLLL_WITH_QD
template<> inline const char* num_type_str<dd_real>()       {return "dd_real";}
template<> inline const char* num_type_str<qd_real>()       {return "qd_real";}
#endif
template<> inline const char* num_type_str<mpfr_t>()      {return "mpfr_t";}


FPLLL_END_NAMESPACE

#endif
