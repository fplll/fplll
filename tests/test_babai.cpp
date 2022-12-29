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
  // Some platforms have sizeof(double) == sizeof(long double)
  // because long double is only required to be at least as large
  // as a double. This means the behaviour of the first test
  // depends on the platform.
  status += test_intrel<mpz_t, long double>(10, 60, sizeof(double) == sizeof(long double));
  status += test_intrel<mpz_t, long double>(10, 70, true);
#endif
#ifdef FPLLL_WITH_QD
  status += test_intrel<mpz_t, dd_real>(10, 110);
  status += test_intrel<mpz_t, dd_real>(10, 120, true);
  // status += test_intrel<mpz_t, qd_real>(10, 100);
  // status += test_intrel<mpz_t, qd_real>(10, 230, true);
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
