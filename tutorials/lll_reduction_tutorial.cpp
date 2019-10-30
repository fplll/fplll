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
 * This function performs LLL reduction on a given lattice base, using the PROVED method.
 * In order to intercat with the LLLReduction class, which implements all LLL reduction methods,
 * the Wrapper interface is used.
 */

void lll_reduction_proved (ZZ_mat<mpz_t> &base) {
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
	cout << base << endl;
}

void lll_reduction_fast(ZZ_mat<mpz_t> &base) {
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
	cout << base << endl;
}

void lll_reduction_heuristic (ZZ_mat<mpz_t> &base) {
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
	cout << base << endl;
}

int main (int argc, char* argv[]) {
	/* We are using a lattice of size 5 * 5,
	 * where the element of each vector is an integer
	 * in the space of [0, 15]. If a different base
	 * is required, the methods for its creation
	 * can be found in  
	 * lattice_base_generation_tutorial.cpp, which
	 * is in the same directory.
	 */
	ZZ_mat<mpz_t> base;
	base.resize(5, 5);
	base.gen_uniform(4);
	cout << base << endl;
	cout << endl;
	lll_reduction_fast(base);
	cout << endl;
	/* Clear and reinitialize */
	base.clear();
	base.resize(5, 5);
	base.gen_uniform(4);
	cout << base << endl;
	cout << endl;
	lll_reduction_fast(base);
	cout << endl;
	/* Clear and reinitialize */
	base.clear();
	base.resize(5, 5);
	base.gen_uniform(4);
	cout << base << endl;
	cout << endl;
	lll_reduction_heuristic(base);
	cout << endl;
}