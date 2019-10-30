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
 * g++ -std=c++11 -march=native -o3 gram-schmidt_operation_tutorial.cpp -lfplll -mpfr -lgmp -o <executable_name_here>
 */

/**
 * NOTE: The MatGSO constructor can receive any of a number of flags as parameters. These flags are:
 * GSO_DEFAULT: Default Gram-Scmidt.
 * GSO_INT_GRAM: The coefficients of Gram matrix are computed in exact arithmetic
 * (of type equal to the one declared at object instantiation). Otherwise, floating
 * point arithmetic is used instead. 
 * GSO_ROW_EXPO: With this method, each row of b is normalized by a power of 2 before doing
 * conversion to floating-point, which hopefully avoids some overflows. This option cannot
 * be enabled if enable_int_gram=true (GSO_INT_GRAM) and works only with FT=double and
 * FT=long double (declared at object instantiation). It is useless and MUST NOT be used
 * for FT=dpe or FT=mpfr_t.
 * GSO_OP_FORCE_LONG: Affects the behaviour of row_addmul(_we). See the documentation of
 * row_addmul.
 */


/**
 * Function that performs Gram-Schmidt orthogonalization on a given lattice base and
 * extracts the Gram-Schmidt orthogonalized base.
 */

void extract_gram_schmidt_base (ZZ_mat<mpz_t> &base, FP_mat<mpfr_t> &gramBase) 
{
	/**
	 * If this is not empty, any operation done on basis is also done on u_matrix. If it
	 * is initialized with the identity matrix, then multiplying the initial basis with
	 * the transformed matrix gives the current basis.
	 */ 
	ZZ_mat<mpz_t> u_matrix; 
	ZZ_mat<mpz_t> u_matrix_inverted_transposed;
	/**
	 * NOTE: Comment retrieved from gso.pyx at fpylll, written by Martin R. Albrecht <martinralbrecht+fpylll@googlemail.com>.
	 * MatGSO provides an interface for performing elementary operations on a basis and computing its
     * Gram matrix and its Gram-Schmidt orthogonalization.  The Gram-Schmidt coefficients are computed
     * on demand. The object keeps track of which coefficients are valid after each row operation.
	 */
	MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> interface (base, u_matrix, u_matrix_inverted_transposed, GSO_INT_GRAM);
	interface.update_gso();
	Matrix<FP_NR<mpfr_t>> r_matrix = interface.get_r_matrix();
	Matrix<FP_NR<mpfr_t>> mu_matrix = interface.get_mu_matrix();
	cout << r_matrix << endl;
	cout << endl;
	cout << mu_matrix << endl;
	cout << endl;
}


int main (int argc, char* argv[])
{
	/* This is the lattice base */
	ZZ_mat<mpz_t> base;
	/* This will store the Gram-Schmidt orthogonalized lattice base */
	FP_mat<mpfr_t> gramBase;
	base.resize(5, 5);
	/** 
	 * For the purposes of this tutorial, a uniform (5*5) matrix
	 * of integers in the space of [0, 15] is generated.
	 * For alternative types of matrices, please refer to
	 * lattice_base_generation_tutorial.cpp.
	 */
	base.gen_uniform(4);
	cout << base << endl;
	cout << endl;
	gramBase.resize(base.get_rows(), base.get_cols());
	gramBase.fill(0.0);
	extract_gram_schmidt_base(base, gramBase);
	cout << gramBase << endl;
}