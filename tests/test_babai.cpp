/* Copyright (C) 2022 Martin R. Albrecht

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

#include "fplll/fplll.h"
#include <vector>

#include <cfloat>  // Needed for precision macros.

using namespace fplll;

template <class ZT, class FT>
int test_intrel(long n, long bits, bool shouldfail = false, int seed = 0)
{
  gmp_randseed_ui(RandGen::get_gmp_state(), seed);
  mpz_t status;
  mpz_init(status);

  ZZ_mat<ZT> B, U, UT;
  B.resize(n, n + 1);
  B.gen_intrel(bits);

  vector<Z_NR<ZT>> t = vector<Z_NR<ZT>>(n + 1);

  for (long i = 0; i < n; i++)
  {
    mpz_urandomb(status, RandGen::get_gmp_state(), 1);
    if (mpz_cmp_si(status, 1))
    {
      t[0].add(t[0], B[i][0]);
    }
  }
  mpz_clear(status);

  vector<Z_NR<ZT>> v = vector<Z_NR<ZT>>(t);

  lll_reduction(B);

  MatGSO<Z_NR<ZT>, FP_NR<FT>> M(B, U, UT, GSO_DEFAULT);
  M.update_gso();
  M.babai(v, 0, n, false);

  vector<Z_NR<ZT>> w = vector<Z_NR<ZT>>(n + 1);

  Z_NR<ZT> tmp;

  for (long i = 0; i < B.get_rows(); i++)
  {
    for (long j = 0; j < B.get_cols(); j++)
    {
      tmp.mul(v[i], B[i][j]);
      w[j].add(w[j], tmp);
    }
  }

  if ((w[0] == t[0]) - shouldfail)
  {
    return 0;
  }
  else
  {
    std::cerr << "n:" << n << ", bits: " << bits << ", shouldfail:" << shouldfail << std::endl;
    std::cerr << "t:" << t << std::endl;
    std::cerr << "w:" << w << std::endl;
    return 1;
  }
}

int main(int argc, char *argv[])
{
  int status = 0;
  RandGen::init_with_seed(0);
  status += test_intrel<mpz_t, double>(10, 20);
  status += test_intrel<mpz_t, double>(10, 30);
  status += test_intrel<mpz_t, double>(10, 40);
  status += test_intrel<mpz_t, double>(10, 50);
  status += test_intrel<mpz_t, double>(10, 60, true);

#ifdef FPLLL_WITH_LONG_DOUBLE
  // long double has various types of platform and compiler specific behaviour.
  // This makes this test case platform specific.
  // These cases are (at least):
  // (1) Some platforms have sizeof(double) == sizeof(long double), because the C standard
  //     only requires long double to be at least as big as double. This means that
  //     in some situations long double may simply alias double (e.g MSVC).
  // (2) Some platforms make long double an alias for IEEE quadruple precision types. This
  //     may mean that long double is actually rather accurate (more than, say, 60 bits of
  //     precision).
  // (3) On x86-64, some compilers may treat long double as the 80-bit extended precision
  //     x87 type. This means that there's more precision than a regular double. In this situation,
  //     the compiler may also make sizeof(long double) == 16, meaning that we cannot detect this
  //     easily based on the size of long double alone.
  //
  // To circumvent these issues, we check how many elements are in the mantissa of long double,
  // using LDBL_MANT_DIG. Specifically, if LDBL_MANT_DIG == DBL_MANT_DIG, then we are in case (1).
  // If they are different, then if LDBL_MANT_DIG < 70 then the second test should also fail.
  status += test_intrel<mpz_t, long double>(10, 60, LDBL_MANT_DIG == DBL_MANT_DIG);
  status += test_intrel<mpz_t, long double>(10, 70, LDBL_MANT_DIG < 70);
#endif
#ifdef FPLLL_WITH_QD
  // QD needs to have the round-to-double flag set on x86 systems. This function
  // does nothing on non-x86 machines (see the readme for the QD library for more).
  // This is also already done in the wrapper.
  unsigned old_cw;
  fpu_fix_start(&old_cw);
  status += test_intrel<mpz_t, dd_real>(10, 110);
  status += test_intrel<mpz_t, dd_real>(10, 120, true);
  // status += test_intrel<mpz_t, qd_real>(10, 100);
  // status += test_intrel<mpz_t, qd_real>(10, 230, true);
  fpu_fix_end(&old_cw);
#endif
  FP_NR<mpfr_t>::set_prec(100);
  status += test_intrel<mpz_t, mpfr_t>(10, 100);
  FP_NR<mpfr_t>::set_prec(200);
  status += test_intrel<mpz_t, mpfr_t>(10, 200);

  if (status == 0)
  {
    std::cerr << "All tests passed." << std::endl;
    return 0;
  }
  else
  {
    return -1;
  }
};
