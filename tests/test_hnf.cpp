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
// #include <fplll.h>
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

template <class ZT> int test_filename(const char *input_filename, HNFMethod method, Z_NR<ZT> det)
{
  int status = 0;
  ZZ_mat<ZT> A, tmp;

  // read the matrix
  status |= read_matrix(A, input_filename);
  if (status)
  {
    cerr << "problem when reading matrix\n";
  }
  tmp = A;

  HNFReduction<ZT> B(A, method, det);

  // one on success
  if (B.is_reduced())
  {
    cerr << "is_hnf_reduced reports success when it should not" << endl;
    return 1;
  }

  B.hnf();

  if (!B.is_reduced())
  {
    cerr << "Output of HNF reduction is not HNF reduced with method " << HNF_METHOD_STR[method]
         // << " test with matrix returned " << get_red_status_str(status) <<
         // endl;
         << " test with matrix returned " << get_red_status_str(B.status) << endl;
  }

  // test if reduction was correct : one on success
  return B.in_lattice_given_hnf(tmp);
}

/**
   @brief Test HNF for a random d x (d+1) matrix with bit size b.

   @param d                dimension
   @param b                bit size
   @param method           HNF method to test

   @return zero on success
*/

template <class ZT> int test_uniform(int d, int b, HNFMethod method)
{
  ZZ_mat<ZT> A, tmp;
  A.resize(d, d + 1);
  A.gen_uniform(b);

  tmp = A;

  Z_NR<ZT> det;
  det = 0;

  HNFReduction<ZT> B(A, method, det);

  B.hnf();
  int status = B.is_reduced();

  if (!status)
  {
    cerr << "Output of HNF reduction is not HNF reduced with method " << HNF_METHOD_STR[method]
         // << " test with matrix returned " << get_red_status_str(status) <<
         // endl;
         << " test with matrix returned " << get_red_status_str(status) << endl;
  }

  // test if reduction was correct : zero on failure
  status = B.in_lattice_given_hnf(tmp);

  if (!status)
  {
    cerr << "Output of HNF reduction is not HNF reduced with method " << HNF_METHOD_STR[method]
         // << " test with matrix returned " << get_red_status_str(status) <<
         // endl;
         << ", basis has changed to a different lattice " << endl;
  }

  return status;
}

/* TODO : KAT tests */

int main(int /*argc*/, char ** /*argv*/)
{
  Z_NR<mpz_t> det;

  int status = 1;
  status &= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_CLASSIC, det);
  status &= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_XGCD, det);
  status &= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_MINORS, det);
  // status &= test_filename(TESTDATADIR "/tests/lattices/dim55_in", HM_PERNETSTEIN, det);

  status &= test_uniform<mpz_t>(15, 10, HM_CLASSIC);
  status &= test_uniform<mpz_t>(15, 10, HM_XGCD);
  status &= test_uniform<mpz_t>(15, 10, HM_MINORS);
  // status &= test_uniform<mpz_t>(15, 10, HM_PERNETSTEIN);

  status &= test_uniform<mpz_t>(20, 5, HM_CLASSIC);
  status &= test_uniform<mpz_t>(20, 5, HM_XGCD);
  status &= test_uniform<mpz_t>(20, 5, HM_MINORS);
  // status &= test_uniform<mpz_t>(20, 5, HM_PERNETSTEIN);

  status &= test_filename(TESTDATADIR "/tests/lattices/example2_in", HM_CLASSIC, det);
  status &= test_filename(TESTDATADIR "/tests/lattices/example2_in", HM_XGCD, det);
  status &= test_filename(TESTDATADIR "/tests/lattices/example2_in", HM_MINORS, det);
  // status &= test_filename(TESTDATADIR "/tests/lattices/example2_in, HM_PERNETSTEIN, det);

  if (status)
  {
    cerr << "All tests passed." << endl;
    return 0;
  }
  else
  {
    return -1;
  }
}
