#ifndef FPLLL_NR_H
#define FPLLL_NR_H

/**
 * loading definitions
 */
#include "../defs.h"
#include <iostream>

#include "fplll/nr/nr_rand.inl"

#include "fplll/nr/nr_Z.inl"
#include "fplll/nr/nr_Z_d.inl"
#include "fplll/nr/nr_Z_l.inl"
#include "fplll/nr/nr_Z_mpz.inl"

#include "fplll/nr/nr_FP.inl"
#include "fplll/nr/nr_FP_d.inl"
#include "fplll/nr/nr_FP_ld.inl"

#ifdef FPLLL_WITH_DPE
#include "fplll/nr/nr_FP_dpe.inl"
#endif

#ifdef FPLLL_WITH_QD
#include "fplll/nr/nr_FP_dd.inl"
#include "fplll/nr/nr_FP_qd.inl"
#endif

#include "fplll/nr/nr_FP_mpfr.inl"

#include "fplll/nr/nr_FP_misc.inl"
#include "fplll/nr/nr_Z_misc.inl"

FPLLL_BEGIN_NAMESPACE

/**
 * Some defaults.
 */

/* Floating-point type inside the SVP/CVP solver */
typedef double enumf;
typedef double enumxt;
// typedef int enumxt;

/**
 * return type
 */
template <class T> inline const char *num_type_str() { return ""; }
#ifdef FPLLL_WITH_ZLONG
template <> inline const char *num_type_str<long>() { return "long"; }
#endif
template <> inline const char *num_type_str<double>() { return "double"; }
template <> inline const char *num_type_str<mpz_t>() { return "mpz_t"; }
#ifdef FPLLL_WITH_LONG_DOUBLE
template <> inline const char *num_type_str<long double>() { return "long double"; }
#endif
#ifdef FPLLL_WITH_DPE
template <> inline const char *num_type_str<dpe_t>() { return "dpe_t"; }
#endif
#ifdef FPLLL_WITH_QD
template <> inline const char *num_type_str<dd_real>() { return "dd_real"; }
template <> inline const char *num_type_str<qd_real>() { return "qd_real"; }
#endif
template <> inline const char *num_type_str<mpfr_t>() { return "mpfr_t"; }

FPLLL_END_NAMESPACE

#endif
