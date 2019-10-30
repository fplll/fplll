#include <stdio.h>
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

void strategize(vector<Strategies> &strategy, block_size)
{
	for (int i = 0; i <= block_size; i++) 
	{
		
	}
}

/**
 * Function that implements BKZ reduction WITHOUT an input pruning vector. Block_size is set to the integer
 * nearest to the square root of the number of vectors in the lattice base. For this, round() function is used.
 */

void bkz_reduction_default (ZZ_mat<mpz_t> &base)
{
	bool status = 0;
	int block_size = round(sqrt(base.get_cols())), precision = 0; /* Block size is set as described above */
	/**
	 * This function implements the BKZ reduction. BKZ_DEFAULT indicates that the default method will be
	 * used (no pruning vector). FT_DOUBLE sets the real number type to double, necessary for Gram-Schmidt step.
	 * Other options are FP_MPFR and FT_DEFAULT. In the case of FT_MPFR, the precision (now set to 0)
	 * MUST be other than 0 (128 is recommended, as it is used on test files).
	 */
	status = bkz_reduction(base, block_size, BKZ_DEFAULT, FT_DOUBLE, precision);
	cout << base << endl;
}
/**
 * Function that implements BKZ reduction WITH an input pruning vector. Block_size is as with bkz_reduction_defalut().
 */

void bkz_reduction_with_parameters (ZZ_mat<mpz_t> &base)
{
	bool status = 0;
	int block_size = round(sqrt(base.get_cols())), precision = 0; /* Block size is set as described in bkz_reduction_default(). */
	/**
	 * Here a vector of Strategies Class is prepared. A strategy covers pruning parameters and
	 * preprocessing block_sizes. This will be our pruning vector. For the purposes of this
	 * tutorial, the same strategy as the one used in test_bkz.cpp will be used here.
	 */
	vector<Strategies> strategy;
	strategize(strategy, block_size);
	status = bkz_reduction(base)
}


 int main (int argc, char *argv[])
 {
 	ZZ_mat<mpz_t> base;
 	base.resize(5, 5);
 	base.gen_uniform(4);
 	cout << base << endl;
 	cout << endl;
 	bkz_reduction_default (base);
 	cout << endl;
 	return 0;
 }