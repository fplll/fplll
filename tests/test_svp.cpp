/* Copyright (C) 2015 Martin Albrecht
 *               2016 Michael Walter

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

enum Test 
{
  SVP_ENUM,
  DSVP_ENUM,
  DSVP_REDUCE
};

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
   @brief Compute the norm of a dual vector (specified by coefficients in the dual basis).

   @param A              input lattice
   @param b              coefficients of shortest dual vector
   @return
*/
template<class ZT>
int dualLength(Float &norm, ZZ_mat<ZT> &A, const IntVect &coords) {
  int d = coords.size();
  if (A.getRows() != d) {
    cerr << "DSVP length error: Coefficient vector has wrong dimension: ";
    cerr << A.getRows() << " vs " << d << endl;
    return 1;
  }
  FloatVect coords_d(d);
  for (int i = 0; i < d; i++) {
    coords_d[i] = coords[i].get_d();
  }
  
  IntMatrix emptyMat;
  MatGSO<Integer, Float> gso(A, emptyMat, emptyMat, GSO_INT_GRAM);
  if (!gso.updateGSO()) {
    cerr << "GSO Failure." << endl;
    return 1;
  }
  Float tmp;
  gso.getR(tmp, d-1, d-1);
  tmp.pow_si(tmp, -1);
  
  FloatVect alpha(d);
  Float mu, alpha2, r_inv;
  norm = 0;
  for (int i = 0; i < d; i++) {
    alpha[i] = coords_d[i];
    for (int j = 0; j < i; j++) {
      gso.getMu(mu, i, j);
      alpha[i] -= mu*alpha[j];
    }
    gso.getR(r_inv, i, i);
    r_inv.pow_si(r_inv, -1);
    alpha2.pow_si(alpha[i], 2);
    norm += alpha2 * r_inv;
  }
  
  return 0;
}

/**
   @brief Test if dual SVP function returns vector with right norm.

   @param A              input lattice
   @param b              shortest dual vector
   @return
*/

template<class ZT>
int testDualSVP(ZZ_mat<ZT> &A, IntVect &b) {
  IntVect solCoord;  // In the LLL-reduced basis
  IntVect solution;
  IntMatrix u;

  Float normb;
  if (dualLength(normb, A, b)) {
    return 1;
  }
  
  int status = lllReduction(A, u, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);
  if (status != RED_SUCCESS) {
    cerr << "LLL reduction failed: " << getRedStatusStr(status) << endl;
    return status;
  }

  status = shortestVector(A, solCoord, SVPM_FAST, SVP_DUAL);

  if (status != RED_SUCCESS) {
    cerr << "Failure: " << getRedStatusStr(status) << endl;
    return status;
  }
  
  Float normSol;
  if (dualLength(normSol, A, solCoord)) {
    return 1;
  }
  
  Float error; error = 1;
  error.mul_2si(error, -(int)error.getprec());
  normb += error;
  if (normSol > normb) {
    cerr << "Returned dual vector too long by more than " << error << endl;
    return 1;
  }

  return 0;
}


/**
   @brief Test if dual SVP reduction returns reduced basis.

   @param A              input lattice
   @param b              shortest dual vector
   @return
*/
template<class ZT>
int testDSVPReduce(ZZ_mat<ZT> &A, IntVect &b) {
  IntMatrix u;
  int d = A.getRows();
  
  Float normb;
  if (dualLength(normb, A, b)) {
    return 1;
  }
  
  int status = lllReduction(A, u, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);
  if (status != RED_SUCCESS) {
    cerr << "LLL reduction failed: " << getRedStatusStr(status) << endl;
    return status;
  }
  
  
  IntMatrix emptyMat;
  MatGSO<Integer, Float> gso(A, emptyMat, emptyMat, GSO_INT_GRAM);
  LLLReduction<Integer, Float> lllObj(gso, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
  
  vector<Strategy> strategies;
  BKZParam dummy(d, strategies);
  BKZReduction<Float> bkzObj(gso, lllObj, dummy);
  bool clean = true;
  
  bkzObj.svp_reduction_ex(0, d, dummy, clean, true);
  status = bkzObj.status;
  if (status != RED_SUCCESS) {
    cerr << "Failure: " << getRedStatusStr(status) << endl;
    return status;
  }
  
  Float normSol;
  Integer zero; zero = 0;
  IntVect e_n(d, zero);
  e_n[d-1] = 1;
  if (dualLength(normSol, A, e_n)) {
    return 1;
  }
  
  Float error; error = 1;
  error.mul_2si(error, -(int)error.getprec());
  normb += error;
  if (normSol > normb) {
    cerr << "Last dual vector too long by more than " << error << endl;
    return 1;
  }

  return 0;
}

/**
   @brief Test if SVP function returns vector with right norm.

   @param input_filename   filename of an input lattice
   @param output_filename  filename of a shortest vector
   @return
*/

template<class ZT>
int testFilename(const char *input_filename, const char *output_filename, const Test test = SVP_ENUM) {
  ZZ_mat<ZT> A;
  readMatrix(A, input_filename);

  IntVect b;
  readVector(b, output_filename);
  
  switch(test) {
    case SVP_ENUM     : return testSVP<ZT>(A, b);
    case DSVP_ENUM    : return testDualSVP<ZT>(A, b);
    case DSVP_REDUCE  : return testDSVPReduce<ZT>(A, b);
  }
  
  cerr << "Unknown test." << endl;
  return 1;
}

/**
   @brief Run SVP tests.

   @param argc             ignored
   @param argv             ignored
   @return
*/

int main(int argc, char *argv[]) {

  int status = 0;
  status |= testFilename<mpz_t>("lattices/example_svp_in", "lattices/example_svp_out");
  status |= testFilename<mpz_t>("lattices/example_dsvp_in", "lattices/example_dsvp_out", DSVP_ENUM);
  status |= testFilename<mpz_t>("lattices/example_dsvp_in", "lattices/example_dsvp_out", DSVP_REDUCE);
  
  if (status == 0) {
    cerr << "All tests passed." << endl;
    return 0;
  } else {
    return -1;
  }

  return 0;
}
