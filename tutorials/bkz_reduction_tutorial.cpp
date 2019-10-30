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

/** 
 * NOTE: We are using BKZ_DEFAULT as the BKZ method. Additional flags (and a combination thereof, with |)
 * that can be used as well are:
 * BKZ_VERBOSE: print additional information during reduction
 * BKZ_NO_LLL: do not run LLL before block reduction (use at your own risk)
 * BKZ_MAX_LOOPS: terminate after max_loops iterations
 * BKZ_MAX_TIME: terminate after max_time time
 * BKZ_BOUNDED_LLL: only run LLL in current block during SVP preprocessing (use at your own risk)
 * BKZ_AUTO_ABORT: heuristically terminate the reduction if progress stalls
 * BKZ_DUMP_GSO: after every iteration write the shape of the current basis to a file
 * BKZ_GH_BND: use the Gaussian heuristic to reduce the enumeration bound of possible
 * BKZ_SD_VARIANT: run SD-BKZ
 * BKZ_SLD_RED: run slide reduction 
 */


void strategize(int block_size, vector<Strategy> &strategies)
{
	for (int i = 0; i <= block_size; i++) 
	{
		Strategy strategy = Strategy::EmptyStrategy(i);
		if ((i == 10) || (i == 20) || (i == 30)) 
		{
			strategy.preprocessing_block_sizes.emplace_back(i / 2);
		}
		strategies.emplace_back(move(strategy));
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
	 * used. FT_DOUBLE sets the real number type to double, necessary for Gram-Schmidt step.
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
	vector<Strategy> strategies;
	strategize(block_size, strategies);
	/* This class generates the BKZ parameters */
	BKZParam parameters (block_size, strategies);
	/**
	 * The flag for the parameter is set to BKZ_DEFAULT. 
	 */
	parameters.flags = BKZ_DEFAULT; 
	status = bkz_reduction(&base, NULL, parameters, FT_DEFAULT, 53);
	cout << base << endl;
}


 int main (int argc, char *argv[])
 {
 	ZZ_mat<mpz_t> base;
 	base.resize(30, 30);
 	base.gen_uniform(10);
 	cout << base << endl;
 	cout << endl;
 	bkz_reduction_default (base);
 	bkz_reduction_with_parameters(base);
 	cout << endl;
 	return 0;
 }