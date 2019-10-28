#include <iostream>
#include <math.h>
#include <mpfr.h>
#include <gmp.h>
#include <fplll.h>
#include <fplll/nr/nr.h>
#include "fplll/defs.h"
#include "fplll/util.h"

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

/* Compile like so: g++ -std=c++11 -march=native -o3 lattice_base_generation_tutorial.cpp -lfplll -lmpfr -lgmp -o <executable_name> 

/** 
 * This function generates a matrix of rows x (rows + 1), similar to 'r' method of latticegen.
 * The first element of each vector is initialized
 * with a random number between 0 and 2^bits - 1. The rest is the canonical vector.
*/
void GenerateKnapsackLatticeBase (ZZ_mat<mpz_t> &base, int rows, bool random, bool with_seed, int random_seed, int bits) {
	/**
	 * This checks whether we desire a random base to be generated.
	 * If so, FPLLL creates a lattice base of random vectors, otherwise it
	 * creates a specific lattice base.
	*/
	int columns = rows;
	if (random) {
		/* If we desire to generate random lattice with a specific seed, we use this function, passing the seed as argument */ 
		if (with_seed){

			RandGen::init_with_seed(random_seed); /* RandGen method initializes a random seed provided by the user */
		}
		else{
			RandGen::init_with_time(); /* RandGen method initializes with current clock time */
		}
		
	}
	/** 
	 * This function resizes the base according to specifications, allocating desired memory.Parameters are:
	 * @param int rows: the number of vectors in matrix.
	 * @param int columns: the number of elements in each vector.
	*/ 
	base.resize(rows, columns + 1);
	/* FPLLL method that gnerates the matrix described above. */
	base.gen_intrel(bits);


}

/**
 * This function generates a matrix of (rows +1) x (rows + 1), similar to 's' method of latticegen.
 * The first vector starts with a random integer in the space of [0, 2^bits-1] and continues with
 * d - 1 independent integers in the space of [0, 2^bits - 1]. For each i > 1, the i-th vector is the
 * i-th canonical unit vector, scaled by a factor of 2^bits.  
*/
void GenerateSimdiophLatticeBase (ZZ_mat<mpz_t> &base, int rows, bool random, bool with_seed, int random_seed, int bits, int bits2) {
	/**
	 * This checks whether we desire a random base to be generated.
	 * If so, FPLLL creates a lattice base of random vectors, otherwise it
	 * creates a specific lattice base.
	 */
	if (random) {
		/* If we desire to generate random lattice base with a specific seed, we use this function, passing the seed as argument */ 
		if (with_seed){

			RandGen::init_with_seed(random_seed); /* RandGen method initializes with random seed */
		}
		else{
			RandGen::init_with_time(); /* RandGen method initializes with current clock time */
		}
		
	}
	/** 
	 * This function resizes the base according to specifications, allocating desired memory.Parameters are:
	 * @param int rows: the number of vectors in matrix.
	 * @param int columns: the number of elements in each vector.
	 */ 
	base.resize(rows +1, rows + 1);
	/* FPLLL method that generates the matrix described above. */
	base.gen_simdioph(bits, bits2);
}

/**
 * Generates a matrix of size dimension, with uniformly random coefficients between 0 and 2^@param int: bit_size - 1
 * The type of randomness used is chosen between clock and pseudorandom of seed equal to @param int: seed
 */

void GenerateUniformlyRandomMatrix (ZZ_mat<mpz_t> &mat, int dimension, int bit_size, bool clock, bool random, int seed){
  if (random) {
      if (clock) { /* If clock randomness is chosen, initializes the random generator to clock */
        RandGen::init_with_time();
      }
      else { /* If pseudorandom randomness is chosen, initializes the random generator to seed with the given seed */
        RandGen::init_with_seed(seed);

      }
  }
  mat.resize(dimension, dimension); /* Resize the matrix to a dimension*dimension matrix */
  mat.gen_uniform(bit_size); /* Generates uniformly random coefficients for the matrix */
}

int main (int argc, char** argv) {
	int seed = 25;
	int dimension = 5;
	int bits2 = 4;
	int bits = 4;
	
	/* This is the lattice Base declaration. */ 
	
	/* ZZ_mat<mpz_t> integer_lattice_base: A matrix of integers. Initially, the matrix is empty. */
	
	ZZ_mat<mpz_t> integer_lattice_base;
	GenerateKnapsackLatticeBase (integer_lattice_base, dimension, false, false, seed, bits);
	cout << "Knapsack-like matrix of dimension " << dimension << " and " << bits << " bits with random seed from clock time" << endl;
	cout << endl;
	cout << integer_lattice_base << endl;
	integer_lattice_base.clear();
	cout << endl;
	GenerateSimdiophLatticeBase (integer_lattice_base, dimension, false, false, seed, bits, bits2);
	cout << "Matrix of form similar to that which is involved in finding rational approximations to reals with the same small denominator, of dimension " << dimension
	<< " with each vector starting with integer of maximum bit-length " << bits2 << " and continues with " << dimension - 1 << " independent integers of maximum bit-length "
	<< bits << ". Each subsequent vector is the canonical unit vector, scaled by a factor of " << pow(bits, 2) << "." << endl;
	cout << endl;
	cout << integer_lattice_base << endl;
	integer_lattice_base.clear();
	cout << endl;
	GenerateUniformlyRandomMatrix (integer_lattice_base, dimension, bits, false, false, seed);
	cout << "Matrix of dimension " << dimension << " whose entries are independent integers of " << bits << " maximum bit-length." << endl;
	cout << endl;
	cout << integer_lattice_base << endl;
	integer_lattice_base.clear();
	cout << endl;
	return 0;
} 
