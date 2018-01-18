/* Copyright (C) 2018 Arnaud Sipasseuth (inspired from the test_lll.cpp)

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
*/

int test_test(ZZ_mat<mpz_t> &A)
{
  // returns 0 if A is a HNF
  int is_reduced = is_hnf_reduced(A);

  if (is_reduced)
    cerr << "is_hnf_reduced reports success when it should not" << endl;

  return is_reduced;
}

/**
   @brief Test HNF reduction.

   @param A                test matrix
   @param method           HNF method to test

   @return zero on success.
*/

int test_hnf(ZZ_mat<mpz_t> &A, HNFMethod method)
{

  ZZ_mat<mpz_t> tmp = A;

  int status = 0;

  // zero on success
  if (test_test(A))
  {
    return 1;
  }

  // zero on success
  status = hnf(A, method);

  // a bit obsolete since there's no error case (yet), but we leave it for the future
  if (status != RED_SUCCESS)
  {
    cerr << "HNF reduction failed with error '" << get_red_status_str(status) << endl;
    return status;
  }

  // zero on success
  status = in_hnf(A, tmp);

  if (status)
    cerr << "Output of HNF reduction is not HNF reduced with method " << HNF_METHOD_STR[method]
         << " test with matrix returned " << status << endl;

  return status;
}

/**
   @brief Test HNF for matrix stored in file pointed to by `input_filename`.

   @param input_filename   a path
   @param method           HNF method to test

   @return zero on success
*/

int test_filename(const char *input_filename, HNFMethod method)
{
  ZZ_mat<mpz_t> A;
  int status = 0;
  status |= read_matrix(A, input_filename);
  status |= test_hnf(A, method);
  return status;
}

/**
   @brief Construct d × (d+1) integer relations matrix with bit size b and test HNF.

   @param d                dimension
   @param b                bit size
   @param method           HNF method to test

   @return zero on success
*/

int test_int_rel(int d, int b, HNFMethod method)
{
  ZZ_mat<mpz_t> A;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  return test_hnf(A, method);
}

int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;
  // status |= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_CLASSIC);
  status |= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_XGCD);
  // status |= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_PERNETSTEIN);

  // status |= test_int_rel(50, 1000, HM_CLASSIC);
  status |= test_int_rel(50, 1000, HM_XGCD);
  // status |= test_int_rel(50, 1000, HM_PERNETSTEIN);

  // status |= test_int_rel(30, 2000, HM_CLASSIC);
  status |= test_int_rel(30, 2000, HM_XGCD);
  // status |= test_int_rel(30, 2000, HM_PERNETSTEIN);

  // status |= test_filename(TESTDATADIR "/tests/lattices/example_out", HM_CLASSIC);
  status |= test_filename(TESTDATADIR "/tests/lattices/example_out", HM_XGCD);
  // status |= test_filename(TESTDATADIR "/tests/lattices/example_out", HM_PERNETSTEIN);

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
