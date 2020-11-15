/* Copyright (C) 2008-2011 Xavier Pujol
   (C) 2015 Michael Walter.
   (C) 2016 Marc Stevens. (generic improvements, auxiliary solutions, subsolutions)
   (C) 2016 Guillaume Bonnoron. (CVP improvements)

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

#include "enumerate.h"
#include <chrono>

FPLLL_BEGIN_NAMESPACE

template class Enumeration<Z_NR<mpz_t>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class Enumeration<Z_NR<mpz_t>, FP_NR<long double>>;
#endif

#ifdef FPLLL_WITH_QD
template class Enumeration<Z_NR<mpz_t>, FP_NR<dd_real>>;

template class Enumeration<Z_NR<mpz_t>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class Enumeration<Z_NR<mpz_t>, FP_NR<dpe_t>>;
#endif

template class Enumeration<Z_NR<mpz_t>, FP_NR<mpfr_t>>;

template class Enumeration<Z_NR<long>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class Enumeration<Z_NR<long>, FP_NR<long double>>;
#endif

#ifdef FPLLL_WITH_QD
template class Enumeration<Z_NR<long>, FP_NR<dd_real>>;

template class Enumeration<Z_NR<long>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class Enumeration<Z_NR<long>, FP_NR<dpe_t>>;
#endif

template class Enumeration<Z_NR<long>, FP_NR<mpfr_t>>;

FPLLL_END_NAMESPACE
