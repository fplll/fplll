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

// Return 0 if all goes well, 1 instead
template <class ZT, class FT> int test_hlll(ZZ_mat<ZT> &A, int flags)
{
  // Define empty u and ut.
  ZZ_mat<ZT> u(0, 0);
  ZZ_mat<ZT> ut(0, 0);
  // Create the MatHouseholder object
  MatHouseholder<Z_NR<ZT>, FP_NR<FT>> Mhouseholder(A, u, ut, flags);
  // Create the HLLLReduction object, on which will be performed the hlll reduction
  HLLLReduction<Z_NR<ZT>, FP_NR<FT>> hlll_obj(Mhouseholder, 0.99, 0.52, 0.99, 0.01, LLL_DEFAULT);

  // The matrix A used to build Mhouseholder is not reduced. If the test return true, then
  // is_hlll_reduced is badly
  // implemented
  int status = is_hlll_reduced<Z_NR<ZT>, FP_NR<FT>>(Mhouseholder, 0.99, 0.52);
  if (status == true)
  {
    cerr << "is_hlll_reduced reports success when it should not." << endl;
    return 1;
  }

  // Perform the hlll reduction
  hlll_obj.lll();

  // Verify if A is hlll reduced thanks to mpfr
  MatHouseholder<Z_NR<ZT>, FP_NR<mpfr_t>> M(A, u, ut, HOUSEHOLDER_DEFAULT);

  // This times, M must be hlll reduced
  status = is_hlll_reduced<Z_NR<ZT>, FP_NR<mpfr_t>>(M, 0.99, 0.52);
  if (status == false)
  {
    cerr << "Output of HLLL reduction is not HLLL reduced." << endl;
    return 1;
  }

  return 0;
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
  // Read A from input_filename
  status |= read_file(A, input_filename);
  // Test hlll reduction on A
  status |= test_hlll<ZT, FT>(A, flags);
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

  int status       = 0;
  int flags_fast   = HOUSEHOLDER_ROW_EXPO /* | HOUSEHOLDER_OP_FORCE_LONG */;
  int flags_proved = HOUSEHOLDER_DEFAULT;

  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/dim55_in", flags_fast);
#ifdef FPLLL_WITH_LONG_DOUBLE
  status |= test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/dim55_in", flags_fast);
#endif  // FPLLL_WITH_LONG_DOUBLE
#ifdef FPLLL_WITH_DPE
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/dim55_in", flags_proved);
#endif  // FPLLL_WITH_DPE
#ifdef FPLLL_WITH_QD
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/dim55_in", flags_fast);
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/dim55_in", flags_fast);
#endif  // FPLLL_WITH_QD
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/dim55_in", flags_proved);

  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                         flags_fast);
#ifdef FPLLL_WITH_LONG_DOUBLE
  status |= test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                              flags_fast);
#endif  // FPLLL_WITH_LONG_DOUBLE
#ifdef FPLLL_WITH_DPE
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                        flags_proved);
#endif  // FPLLL_WITH_DPE
#ifdef FPLLL_WITH_QD
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                          flags_fast);
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                          flags_fast);
#endif  // FPLLL_WITH_QD
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                         flags_proved);

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
