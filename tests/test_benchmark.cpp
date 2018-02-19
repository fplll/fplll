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
#include <stdlib.h>     /* srand, rand */
#include <gso.h>
#include <gso_gram.h>
#include <gso_interface.h>
#include <gso_givens.h>
#include <nr/matrix.h>
#include <test_utils.h>
//#include <random>

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif


template<class ZT, class FT> int benchmark(Matrix<ZT> A)
{ 

  Matrix<ZT> B,U,UT;
  int rowsA = A.get_rows();
  int colsA = A.get_cols();

  B.resize(rowsA,colsA);


  for(int i = 0; i < rowsA; i++) {
    for(int j = 0; j < colsA; j++) {
        B(i,j) = A(i,j);
    }
  }

  MatGSO<ZT, FT> M(A, U, UT, GSO_ROW_EXPO);
  MatGSOGivens<ZT, FT> M_givens(B, U, UT,  GSO_ROW_EXPO | GSO_GIVENS_MOVE_LAZY );
  LLLReduction<ZT, FT> LLLObj(M, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_VERBOSE);  
  LLLReduction<ZT, FT> LLLObj_givens(M_givens, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_VERBOSE);

  cerr << "Givens" << endl;
  LLLObj_givens.lll();

  cerr << "GSO" << endl;  
  LLLObj.lll();  

  M.update_gso();
  M_givens.update_gso();


  int is_reduced = is_lll_reduced<ZT, FT>(M, LLL_DEF_DELTA, LLL_DEF_ETA);

  int is_reduced_givens = is_lll_reduced<ZT, FT>(M_givens, LLL_DEF_DELTA, LLL_DEF_ETA);

    if (is_reduced != 1)
    {
      cerr << "::::::::::::::The basis GSO-object is not LLL-reduced after calling LLL\n";
      return 1;
    }
    if (is_reduced_givens != 1)
    {
      cerr << ":::::::::::::::The givens GSO-object is not LLL-reduced after calling LLL\n";
      return 1;
    }

return 0;


}


template <class ZT, class FT> int test_int_rel(int d, int b)
{
  ZZ_mat<ZT> A;
  A.resize(d,d+1);
  A.gen_intrel(b);

  return benchmark<Z_NR<ZT>, FT>(A);
}



int main(int /*argc*/, char ** /*argv*/)
{
  ZZ_mat<mpz_t> A;
  int status = 0;

  int sequence_length = 8;
  int dimension_sequence[8] = {8, 8, 40, 100, 128, 140, 160, 200};
  int bitsize_sequence[8] = {2000, 200000, 1600, 10000, 10000, 10000, 16000, 20000};
  for(int i = 0; i < sequence_length; i++)
    test_int_rel<mpz_t, FP_NR<double>>(dimension_sequence[i], bitsize_sequence[i]);

  if (status == 0)
  {
    cerr << "All tests passed." << endl;
    return 0;
  }
  else
  {
    return -1;
  }

  return 0;
}
