/* Copyright (C) 2015 Martin Albrecht

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

#include <cstring>
#include <hlll.h>
#include <householder.h>
#include <test_utils.h>

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

template <class ZT, class FT> int test_lll(ZZ_mat<ZT> &A, int flags, int prec = 53)
{
  ZZ_mat<ZT> u(0, 0);
  MatHouseholder<Z_NR<ZT>, FP_NR<FT>> Mhouseholder(A, u, flags);
  HLLLReduction<Z_NR<ZT>, FP_NR<FT>> hlll_obj(Mhouseholder, 0.99, 0.52, 0.99, 0.01, LLL_DEFAULT);

  hlll_obj.lll();

  MatHouseholder<Z_NR<ZT>, FP_NR<mpfr_t>> M(A, u, HOUSEHOLDER_DEFAULT);
  int status = is_hlll_reduced<Z_NR<ZT>, FP_NR<mpfr_t>>(M, 0.99, 0.52);

  if (status == 0)
    cerr << "Output of HLLL reduction is not HLLL reduced." << endl;

  return (status == 0);
}

/**
   @brief Test HLLL for matrix stored in file pointed to by `input_filename`.

   @param input_filename   a path
   @param flags            flags to use

   @return zero on success
*/

template <class ZT, class FT> int test_filename(const char *input_filename, int flags)
{
  ZZ_mat<ZT> A;
  int status = 0;
  status |= read_matrix(A, input_filename);
  status |= test_lll<ZT, FT>(A, flags);
  return status;
}

/**
   @brief Construct d × (d+1) integer relations matrix with bit size b and test LLL.

   @param d                dimension
   @param b                bit size
   @param method           LLL method to test
   @param float_type       floating point type to test
   @param flags            flags to use
   @param prec             precision used for is_lll_reduced

   @return zero on success
*/

int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;

  status |=
      test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/dim55_in", HOUSEHOLDER_ROW_EXPO);
#ifdef FPLLL_WITH_LONG_DOUBLE
  status |= test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/dim55_in",
                                              HOUSEHOLDER_ROW_EXPO);
#endif  // FPLLL_WITH_LONG_DOUBLE
#ifdef FPLLL_WITH_DPE
  status |=
      test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/dim55_in", HOUSEHOLDER_DEFAULT);
#endif  // FPLLL_WITH_DPE
#ifdef FPLLL_WITH_QD
  status |=
      test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/dim55_in", HOUSEHOLDER_ROW_EXPO);
  status |=
      test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/dim55_in", HOUSEHOLDER_ROW_EXPO);
#endif  // FPLLL_WITH_QD
  status |=
      test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/dim55_in", HOUSEHOLDER_DEFAULT);

  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                         HOUSEHOLDER_ROW_EXPO);
#ifdef FPLLL_WITH_LONG_DOUBLE
  status |= test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                              HOUSEHOLDER_ROW_EXPO);
#endif  // FPLLL_WITH_LONG_DOUBLE
#ifdef FPLLL_WITH_DPE
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                        HOUSEHOLDER_DEFAULT);
#endif  // FPLLL_WITH_DPE
#ifdef FPLLL_WITH_QD
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                          HOUSEHOLDER_ROW_EXPO);
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                          HOUSEHOLDER_ROW_EXPO);
#endif  // FPLLL_WITH_QD
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                         HOUSEHOLDER_DEFAULT);

  if (status == 0)
  {
    cerr << "All tests passed." << endl;
    return 0;
  }
  else
  {
    return -1;
  }
}
