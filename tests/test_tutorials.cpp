/*

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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <fplll.h>
#include <fplll/nr/nr.h>
#include <fplll/wrapper.h>
#include "fplll/defs.h"
#include "fplll/util.h"

using namespace std;
using namespace fplll;

/**
   @brief Read T from `input_filename`.
   @param X T (T is usually a ZZ_mat<ZT> or a vector<Z_NR<ZT>>
   @param input_filename
   @return zero if the file is correctly read, 1 otherwise.
*/

template <class T> int read_file(T &X, const char *input_filename) {
  int status = 0;
  ifstream is;
  is.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    is.open(input_filename);
    is >> X;
    is.close();
  }
  catch (const ifstream::failure&) {
    status = 1;
    cerr << "Error by reading " << input_filename << "." << endl;
    cout << is.rdstate() << endl;
  }

  return status;
}

/* @brief test base generation tutorial
   @param ZZ_mat<ZT> lattice_base:   input lattice base.
   @param int rows:                  number of vectors in base.
   @param int columns:               number of elements in each vector.
   @param string type:               type of generation to be tested;
   
   @return 0 on success, where the input matrix will be compared with pre-initialized matrix according to method.
*/
template <class ZT> bool test_generation_tutorial (ZZ_mat<ZT> &lattice_base, const char *input_name) 
{
   ZZ_mat<ZT> test_matrix;
   int status = 1, input = 1, rows = lattice_base.get_rows(), columns = lattice_base.get_cols();
   input = read_file(test_matrix, input_name);
   if (input != 1)
   {
      for (int i = 0; i < rows; i++)
      {
         for (int j = 0; j < columns; j++)
         {
            status = test_matrix[i][j] != lattice_base[i][j];
         }
      }
   }
   return status;

}

/* @brief Function that clears and resizes base.
   @param ZZ_mat<ZT> input_base: input lattice base.
   @param int rows:              number of vectors in base.
   @param int columns:           number of elements in each vector.

*/

template <class ZT> void clear_base (ZZ_mat<ZT> &input_base, int rows, int columns) 
{
   input_base.clear();
   input_base.resize(rows, columns);
}

/* @brief Function that clears and resizes base.
   @param ZZ_mat<ZT> input_base: input lattice base.
   @param int rows:              number of vectors in base.
   @param int columns:           number of elements in each vector.

*/

template <class ZT> int test_generation_tutorial (ZZ_mat<ZT> &test_base)
{
   int status = 0;
   test_base.resize(5, 6);
   test_base.gen_intrel(4);
   status |= test_generation_tutorial (test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/intrel_base_example");
   clear_base(test_base, 6, 6);
   test_base.gen_simdioph(4, 4);
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/simdioph_base_example");
   clear_base(test_base, 5, 5);
   test_base.gen_uniform(4);
   status  |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/uniform_base_example");
   clear_base(test_base, 10, 10);
   test_base.gen_ntrulike(4);
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/ntrulike_base_example");
   clear_base(test_base, 10, 10);
   test_base.gen_ntrulike_withq(4);
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/ntrulike_withq_base_example");
   clear_base(test_base, 10, 10);
   test_base.gen_ntrulike2(4);
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/ntrulike2_base_example");
   clear_base(test_base, 10, 10);
   test_base.gen_ntrulike2_withq(4);
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/ntrulike2_withq_base_example");
   clear_base(test_base, 5, 5);
   test_base.gen_qary_prime(1, 4);
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/qary_prime_base_example");
   clear_base(test_base, 5, 5);
   test_base.gen_qary_withq(1, 4);
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/qary_withq_base_example");
   clear_base(test_base, 5, 5);
   test_base.gen_qary(1, 4);
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/qary_base_example");
   clear_base(test_base, 5, 5);
   test_base.gen_trg(1.02);
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/trg_base_example");
   FP_NR<mpfr_t> *input_vector = new FP_NR<mpfr_t>[test_base.get_cols()];
   for (int i =0; i < test_base.get_cols(); i++)
   {
      input_vector[i] = (float) i + 1;
   }
   clear_base(test_base, 5, 5);
   test_base.gen_trg2(input_vector);
   delete[] input_vector;
   status |= test_generation_tutorial(test_base, "/home/marios/Desktop/fplll_days_5/fplll/tests/lattices/trg2_base_example");
   return status;
}

int main (int argc, char * argv[])
{
   int status = 0;
   ZZ_mat<mpz_t> test_base;
   status |=test_generation_tutorial(test_base);
   cout << status << endl;
}