/* Copyright (c) 2019 Marios Mavropoulos Papoudas

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
#include <math.h>
#include <mpfr.h>
#include <gmp.h>
#include <fplll.h>
#include <fplll/nr/nr.h>
#include <fplll/wrapper.h>
#include "fplll/defs.h"
#include "fplll/util.h"

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

/**
 * Compile like so:
 * g++ -std=c++11 -march=native -o3 lll_reduction_tutorial.cpp -lfplll -mpfr -lgmp -o <executable_name_here>
 */

/**
 * This function performs LLL reduction on a given lattice base, using the PROVED method.
 * In order to intercat with the LLLReduction class, which implements all LLL reduction methods,
 * the Wrapper interface is used.
 */

void lll_reduction_proved (ZZ_mat<mpz_t> &base)
 {
	bool result = 0;
	double delta = 0.99, eta = 0.51; /* These are the default values of delta and eta respectively. */
	ZZ_mat<mpz_t> identity_matrix; /* This must be either empty or the identity matrix */
	ZZ_mat<mpz_t> identity_matrix_transposed; /* As with identity_matrix */
	/** 
	 * LM_PROVED denotes that the proved method of lll will be used. FT_MPFR denotes the data type
	 * of real numbers, necessary for the Gram-Schmidt step. For our purposes, this is set to FT_MPFR.
	 * Other options are FT_DOUBLE, which sets to double precision, and FT_DEFAULT, which sets to 
	 * default precision.
	 */
	result = lll_reduction(base, identity_matrix, identity_matrix_transposed, delta, eta, LM_PROVED, FT_MPFR, 0, LLL_DEFAULT);
}


/**
   @brief Write T to `output_filename`.
   @param X T (T is usually a ZZ_mat<ZT> or a vector<Z_NR<ZT>>
   @param output_filename
   @return zero if the file is correctly written to, 1 otherwise.
*/

template <class T> int write_to_file(T &X, const char *output_filename) {
  int status = 0;
  ofstream os;
  os.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  try {
    os.open(output_filename, std::ios::app);
    os << X;
    os.close();
  }
  catch (const ofstream::failure&) {
    status = 1;
    cerr << "Error by writing to " << output_filename << "." << endl;
    cout << os.rdstate() << endl;
  }

  return status;
}

void lll_reduction_fast(ZZ_mat<mpz_t> &base)
 {
	bool result = 0;
	double delta = 0.99, eta = 0.51; /* These are the default values of delta and eta respectively. */
	ZZ_mat<mpz_t> identity_matrix; /* This must be either empty or the identity matrix */
	ZZ_mat<mpz_t> identity_matrix_transposed; /* As with identity_matrix */
	/**
	 * LM_FAST denotes fast lll method will be used. The real numbers data type HAS TO BE either double, long double,
	 * double-double or quad-double for the fast method. In this case, it has been set to FT_DOUBLE.
	 * Other options are FT_LONG_DOUBLE, FT_DD and FT_QD. For the last two options, compiler support for Quad Double
	 * and vouble double data types is required. For all other flags, the same rules as with lll_reduction_proved() apply. 
	 */
	result = lll_reduction(base, identity_matrix, identity_matrix_transposed, delta, eta, LM_FAST, FT_DOUBLE, 0, LLL_DEFAULT);
}

void lll_reduction_heuristic (ZZ_mat<mpz_t> &base) 
{
	bool result = 0;
	double delta = 0.99, eta = 0.51; /* These are the default values of delta and eta respectively. */
	ZZ_mat<mpz_t> identity_matrix; /* This must be either empty or the identity matrix */
	ZZ_mat<mpz_t> identity_matrix_transposed; /* As with identity_matrix */
	/**
	 * LM_HEURISTIC denotes that the heuristic lll method will be used. The real numbers data type
	 * HAS TO BE either default precision or DPE for the heuristic method. In this case, it has
	 * been set to FT_DEFAULT. Other options are FT_DPE. For all other flags, the same rules
	 * as with lll_reduction_proved() apply. 
	 */
	result = lll_reduction(base, identity_matrix, identity_matrix_transposed, delta, eta, LM_HEURISTIC, FT_DEFAULT, 0, LLL_DEFAULT);
}

int main (int argc, char* argv[]) 
{
	/* We are using a lattice of size 5 * 5,
	 * where the element of each vector is an integer
	 * in the space of [0, 15]. If a different base
	 * is required, the methods for its creation
	 * can be found in  
	 * lattice_base_generation_tutorial.cpp, which
	 * is in the same directory.
	 */
	ZZ_mat<mpz_t> base;
	int output = 0;
	base.resize(5, 5);
	base.gen_uniform(4);
	lll_reduction_fast(base);
	output = write_to_file(base, "lll_output");
	/* Clear and reinitialize */
	base.clear();
	base.resize(5, 5);
	base.gen_uniform(4);
	lll_reduction_proved(base);
	output = write_to_file(base, "lll_output");
	/* Clear and reinitialize */
	base.clear();
	base.resize(5, 5);
	base.gen_uniform(4);
	lll_reduction_heuristic(base);
	output = write_to_file(base, "lll_output");
	return 0;
}