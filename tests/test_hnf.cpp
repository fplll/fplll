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

/* this is to make tests run... */
#include <../fplll/hnf.cpp>

#include <cstring>
#include <fplll.h>
#include <test_utils.h>

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

/**
   @brief Test HNF for matrix stored in file pointed to by `input_filename`.

   @param input_filename   a path
   @param method           HNF method to test

   @return zero on success
*/

int test_filename(const char *input_filename, HNFMethod method)
{
  ZZ_mat<mpz_t> A, tmp;

  // have the status
  int status = 0;

  // read the matrix
  status |= read_matrix(A, input_filename);
  if (status)
  {
    cerr << "problem when reading matrix\n";
  }
  tmp = A;

  // zero on success
  if (!is_hnf_reduced(A))
  {
    cerr << "is_hnf_reduced reports success when it should not" << endl;
    return 1;
  }

  status = hnf(A, method);

  if (status != RED_SUCCESS)
  {
    cerr << "Output of HNF reduction is not HNF reduced with method " << HNF_METHOD_STR[method]
         // << " test with matrix returned " << get_red_status_str(status) << endl;
         << " test with matrix returned " << get_red_status_str(status) << endl;
  }

  // test if reduction was correct : zero on success
  status = in_hnf(A, tmp);

  return status;
}

int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;
  // status |= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_CLASSIC);
  status |= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_XGCD);
  // status |= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_PERNETSTEIN);

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
