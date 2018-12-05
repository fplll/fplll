/* Copyright (C) 2015 Martin Albrecht
   Copyright (C) 2017-2018 Laurent Grémy

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
#include <fplll.h>
#include <test_utils.h>

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

/**
   @brief Test the tester.

   @param A
   @return zero on success.

   Should return 0, since the at this point of the program, A is not reduced.
*/

template <class ZT> int test_test(ZZ_mat<ZT> &A)
{
  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  MatHouseholder<Z_NR<ZT>, FP_NR<mpfr_t>> M(A, U, UT, 0);

  int is_reduced =
      is_hlll_reduced<Z_NR<ZT>, FP_NR<mpfr_t>>(M, LLL_DEF_DELTA, LLL_DEF_ETA, HLLL_DEF_THETA);

  if (is_reduced == RED_SUCCESS)
    cerr << "is_hlll_reduced reports success when it should not" << endl;

  return (is_reduced == RED_SUCCESS);
}

/**
   @brief Test HLLL reduction.

   @param A                test matrix
   @param method           LLL method to test
   @param float_type       floating point type to test
   @param flags            flags to use

   @return zero on success.
*/

template <class ZT>
int test_hlll(ZZ_mat<ZT> &A, LLLMethod method, FloatType float_type, int flags = LLL_DEFAULT)
{

  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  int status = 0;

  // zero on success
  if (test_test(A))
  {
    return 1;
  }

  // zero on success
  status = hlll_reduction(A, LLL_DEF_DELTA, LLL_DEF_ETA, HLLL_DEF_THETA, HLLL_DEF_C, method,
                          float_type, 0, flags, false);

  if (status != RED_SUCCESS)
  {
    cerr << "HLLL reduction failed with error '" << get_red_status_str(status);
    cerr << "' for method " << LLL_METHOD_STR[method];
    cerr << " and float type " << FLOAT_TYPE_STR[float_type] << endl;
    return status;
  }

  int prec = hlll_min_prec(A.get_rows(), A.get_cols(), LLL_DEF_DELTA, LLL_DEF_ETA, HLLL_DEF_THETA,
                           HLLL_DEF_C);

  const int old_prec = FP_NR<mpfr_t>::set_prec(prec);

  MatHouseholder<Z_NR<ZT>, FP_NR<mpfr_t>> M(A, U, UT, 0);

  // one on success
  status = is_hlll_reduced<Z_NR<ZT>, FP_NR<mpfr_t>>(M, LLL_DEF_DELTA, LLL_DEF_ETA, HLLL_DEF_THETA);

  FP_NR<mpfr_t>::set_prec(old_prec);

  if (status != RED_SUCCESS)
    cerr << "Output of HLLL reduction is not HLLL reduced with method " << LLL_METHOD_STR[method]
         << endl;

  return (status != RED_SUCCESS);
}

/**
   @brief Test HLLL for matrix stored in file pointed to by `input_filename`.

   @param input_filename   a path
   @param method           LLL method to test
   @param float_type       floating point type to test
   @param flags            flags to use

   @return zero on success
*/

template <class ZT>
int test_filename(const char *input_filename, LLLMethod method, FloatType float_type = FT_DEFAULT,
                  int flags = LLL_DEFAULT)
{
  ZZ_mat<ZT> A;
  int status = 0;
  status |= read_file(A, input_filename);
  status |= test_hlll<ZT>(A, method, float_type, flags);
  return status;
}

/**
   @brief Construct d × (d+1) integer relations matrix with bit size b and test LLL.

   @param d                dimension
   @param b                bit size
   @param method           LLL method to test
   @param float_type       floating point type to test
   @param flags            flags to use

   @return zero on success
*/

template <class ZT>
int test_int_rel(int d, int b, LLLMethod method, FloatType float_type = FT_DEFAULT,
                 int flags = LLL_DEFAULT)
{
  ZZ_mat<ZT> A;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  return test_hlll<ZT>(A, method, float_type, flags);
}

int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/dim55_in", LM_WRAPPER, FT_DEFAULT);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/dim55_in", LM_PROVED, FT_MPFR);

  status |= test_int_rel<mpz_t>(50, 1000, LM_FAST, FT_DOUBLE);
  status |= test_int_rel<mpz_t>(30, 1000, LM_PROVED, FT_MPFR);

  status |= test_int_rel<mpz_t>(30, 2000, LM_PROVED, FT_DPE);
  status |= test_int_rel<mpz_t>(20, 2000, LM_PROVED, FT_MPFR);

  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", LM_FAST, FT_DOUBLE);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", LM_PROVED, FT_MPFR);

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
