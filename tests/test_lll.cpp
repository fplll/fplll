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
#include <fplll.h>

using namespace std;
using namespace fplll;

template<class ZT>
void readMatrix(ZZ_mat<ZT> &A, const char *input_filename) {
  istream *is = new ifstream(input_filename);
  *is >> A;
  delete is;
}

/**
   @brief Test the tester.

   @param A
   @return zero on success.
*/

template<class ZT>
int testTest(ZZ_mat<ZT> &A) {
  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  MatGSO<Z_NR<ZT>, FP_NR<mpfr_t> > M(A, U, UT, 0);

  int is_reduced = isLLLReduced<Z_NR<ZT>, FP_NR<mpfr_t> >(M, LLL_DEF_DELTA, LLL_DEF_ETA);

  if (is_reduced)
    cerr << "isLLLReduced reports success when it should not" << endl;

  return is_reduced;
}

/**
   @brief Test LLL reduction.

   @param A                test matrix
   @param method           LLL method to test
   @param float_type       floating point type to test
   @param flags            flags to use

   @return zero on success.
*/

template<class ZT>
int testLLL(ZZ_mat<ZT> &A, LLLMethod method, FloatType float_type, int flags=LLL_DEFAULT) {

  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  int status = 0;

  // zero on success
  if (testTest(A)){
    return 0;

  }

  // zero on success
  status = lllReduction(A, LLL_DEF_DELTA, LLL_DEF_ETA, method, float_type, 0, flags);
  if (status != RED_SUCCESS) {
    cerr << "LLL reduction failed with error '" << getRedStatusStr(status);
    cerr << "' for method " << LLL_METHOD_STR[method];
    cerr << " and float type " << FLOAT_TYPE_STR[float_type] << endl;
    return status;
  }

  MatGSO<Z_NR<ZT>, FP_NR<mpfr_t> > M(A, U, UT, 0);

  // one on success
  status = isLLLReduced<Z_NR<ZT>, FP_NR<mpfr_t> >(M, LLL_DEF_DELTA, LLL_DEF_ETA);

  if (status == 0)
    cerr << "Output of LLL reduction is not LLL reduced with method " << LLL_METHOD_STR[method] << endl;

  return (status == 0);
}

/**
   @brief Test LLL for matrix stored in file pointed to by `input_filename`.

   @param input_filename   a path
   @param method           LLL method to test
   @param float_type       floating point type to test
   @param flags            flags to use

   @return zero on success
*/

template<class ZT>
int testFilename(const char *input_filename, LLLMethod method, FloatType float_type=FT_DEFAULT, int flags=LLL_DEFAULT) {
  ZZ_mat<ZT> A;
  readMatrix(A, input_filename);
  return testLLL<ZT>(A, method, float_type, flags);
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

template<class ZT>
int testIntRel(int d, int b, LLLMethod method, FloatType float_type=FT_DEFAULT, int flags=LLL_DEFAULT) {
  ZZ_mat<ZT> A;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  return testLLL<ZT>(A, method, float_type, flags);
}

int main(int argc, char *argv[]) {

  int status = 0;
  status |= testFilename<mpz_t>("lattices/dim55_in", LM_WRAPPER, FT_DEFAULT);
  status |= testFilename<mpz_t>("lattices/dim55_in", LM_PROVED, FT_MPFR);

  status |= testIntRel<mpz_t>(50, 1000, LM_FAST, FT_DOUBLE);
  status |= testIntRel<mpz_t>(50, 1000, LM_PROVED, FT_MPFR);

  status |= testIntRel<mpz_t>(30, 2000, LM_HEURISTIC, FT_DPE);
  status |= testIntRel<mpz_t>(30, 2000, LM_PROVED, FT_DPE);
  status |= testIntRel<mpz_t>(30, 2000, LM_PROVED, FT_MPFR);


  status |= testFilename<mpz_t>("lattices/example_in", LM_HEURISTIC);
  status |= testFilename<mpz_t>("lattices/example_in", LM_FAST, FT_DOUBLE);
  status |= testFilename<mpz_t>("lattices/example_in", LM_PROVED, FT_MPFR);
  status |= testFilename<mpz_t>("lattices/example_in", LM_FAST, FT_DOUBLE, LLL_DEFAULT | LLL_EARLY_RED);
  status |= testFilename<mpz_t>("lattices/example_in", LM_HEURISTIC, FT_DEFAULT, LLL_DEFAULT | LLL_EARLY_RED);

  if (status == 0) {
    cerr << "All tests passed." << endl;
    return 0;
  } else {
    return -1;
  }

  return 0;
}
