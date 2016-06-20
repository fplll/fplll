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

/**
   @brief Read matrix from `input_filename`.

   @param A
   @param input_filename
   @return
*/

template<class ZT>
void readMatrix(ZZ_mat<ZT> &A, const char *input_filename) {
  istream *is = new ifstream(input_filename);
  *is >> A;
  delete is;
}

/**
   @brief Read vector from `input_filename` into `b`.

   @param b                vector
   @param input_filename   filename
   @return
*/

template<class ZT>
void readVector(vector<Z_NR<ZT> > &b, const char *input_filename) {
  istream *is = new ifstream(input_filename);
  *is >> b;
  delete is;
}

/**
   @brief Test if SVP function returns vector with right norm.

   @param A              input lattice
   @param b              shortest vector
   @return
*/

template<class ZT>
int testSVP(ZZ_mat<ZT> &A, IntVect &b) {
  IntVect solCoord;  // In the LLL-reduced basis
  IntVect solCoord2; // In the initial basis
  IntVect solution;
  IntMatrix u;

  int status = lllReduction(A, u, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);
  if (status != RED_SUCCESS) {
    cerr << "LLL reduction failed: " << getRedStatusStr(status) << endl;
    return status;
  }

  status = shortestVector(A, solCoord, SVPM_PROVED, SVP_DEFAULT);

  if (status != RED_SUCCESS) {
    cerr << "Failure: " << getRedStatusStr(status) << endl;
    return status;
  }

  vectMatrixProduct(solCoord2, solCoord, u);
  vectMatrixProduct(solution, solCoord, A);

  Z_NR<ZT> tmp;
  Z_NR<ZT> norm_s;
  Z_NR<ZT> norm_b;

  for(int i=0; i<A.getCols(); i++) {
    tmp.mul(solution[i], solution[i]);
    norm_s.add(norm_s, tmp);

    tmp.mul(b[i], b[i]);
    norm_b.add(norm_b, tmp);
  }
  if (norm_s != norm_b)
    return 1;

  return 0;
}

/**
   @brief Test if SVP function returns vector with right norm.

   @param input_filename   filename of an input lattice
   @param output_filename  filename of a shortest vector
   @return
*/

template<class ZT>
int testFilename(const char *input_filename, const char *output_filename) {
  ZZ_mat<ZT> A;
  readMatrix(A, input_filename);

  IntVect b;
  readVector(b, output_filename);

  return testSVP<ZT>(A, b);
}

/**
   @brief Run SVP tests.

   @param argc             ignored
   @param argv             ignored
   @return
*/

int main(int argc, char *argv[]) {

  int status = 0;
  if (argc >= 2)
   status |= testFilename<mpz_t>(argv[1] /*"lattices/example_svp_in"*/, "lattices/example_svp_out");
  else
   status |= testFilename<mpz_t>("lattices/example_svp_in", "lattices/example_svp_out");

  if (status == 0) {
    cerr << "All tests passed." << endl;
    return 0;
  } else {
    return -1;
  }

  return 0;
}
