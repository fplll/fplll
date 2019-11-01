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
#include "fplll/defs.h"
#include "fplll/util.h"

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

/* Compile like so: g++ -std=c++11 -march=native -o3 lattice_base_generation_tutorial.cpp -lfplll -lmpfr -lgmp -o <executable_name> */


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

/**
 * This function decides the randomness seed. By default, the same seed is used every time to ensure reproducibility.
 * However, there exists the option of providing a random integer seed or of randomizing using seed provided by CPU clock.
 * This is decided here.
 */

void decide_randomness (bool generate_with_clock) 
{
	int seed = 0;
	if (generate_with_clock) 
	{
		RandGen::init_with_time();
	}
	else 
	{
		cin >> seed;
		RandGen::init_with_seed(seed);
	}
}

/** 
 * This function generates a matrix of rows x (rows + 1), similar to 'r' method of latticegen.
 * The first element of each vector is initialized with a random number between
 * 0 and 2^bit_size - 1. The rest is the canonical vector.
 */

void generate_knapsack_matrix (ZZ_mat<mpz_t> &base, int rows, int bit_size) 
{
	
	int columns = rows;
	/** 
	 * This function resizes the base according to specifications, allocating desired memory.Parameters are:
	 * @param int rows: the number of vectors in matrix.
	 * @param int columns: the number of elements in each vector.
	 */ 
	base.resize(rows, columns + 1);
	/* FPLLL method that gnerates the matrix described above. */
	base.gen_intrel(bit_size);


}

/**
 * This function generates a matrix of (rows +1) x (rows + 1), similar to 's' method of latticegen.
 * The first vector starts with a random integer in the space of [0, 2^bit_size-1] and continues with
 * d - 1 independent integers in the space of [0, 2^bit_size - 1]. For each i > 1, the i-th vector is the
 * i-th canonical unit vector, scaled by a factor of 2^bit_size.  
 */
void generate_simdioph_matrix (ZZ_mat<mpz_t> &base, int rows, int bit_size, int bit_size2) 
{
	/** 
	 * This function resizes the base according to specifications, allocating desired memory.Parameters are:
	 * @param int rows: the number of vectors in matrix.
	 * @param int columns: the number of elements in each vector.
	 */ 
	base.resize(rows +1, rows + 1);
	/* FPLLL method that generates the matrix described above. */
	base.gen_simdioph(bit_size, bit_size2);
}

/**
 * Generates a matrix of size dimension, with uniformly random coefficients between 0 and 2^bit_size - 1
 */
void generate_uniformly_random_matrix (ZZ_mat<mpz_t> &base, int dimension, int bit_size)
{
  base.resize(dimension, dimension); /* Resize the matrix to a dimension * dimension matrix */
  base.gen_uniform(bit_size); /* Generates uniformly random coefficients for the matrix */
}


/**
 * Generates a matrix of size (2 * dimension) * (2 * dimension), as the 'n' method of latticegen API.
 * The method samples a uniform h in the ring Z_q[x]/(x^n-1). It finally returns the 2 x 2 block matrix [[I, Rot(h)], [0, q*I]],
 * where each block is d x d, the first row of Rot(h) is the coefficient vector of h, and the i-th row of Rot(h)
 * is the shift of the (i-1)-th (with last entry put back in first position), for all i>1.
 * Warning: this does not produce a genuine ntru lattice with h a genuine public key.
 */
void generate_ntru_like_matrix (ZZ_mat<mpz_t> &base, int dimension, int bit_size, char selection_char)
{
  base.resize(2 * dimension, 2 * dimension); /* Resize the matrix to a (2 * dimension) * (2 * dimension) matrix */
  /** 
   * If the character is  'b', then the method first samples an integer q in space [0, 2^b-1] 
   */
  if (selection_char == 'b')
  {
  	base.gen_ntrulike(bit_size); 
  }
  else 
   {
  	if (selection_char == 'q') 
  	{
  		base.gen_ntrulike_withq(bit_size); /* If the character is 'q', then the integer is set to the provided value */
  	}
  	else
  	{
  		cout << "Warning: selection character MUST be 'b' or 'q'." << endl;
  	}
  }
}

/**
 * Generates a matrix just like generate_ntrulike_matrix(), except that the constructed matrix is [[q*I, 0], [Rot(h), I]].
 * This functions just like 'N' method of latticegen API.
 */
void generate_ntru_like_matrix_alt (ZZ_mat<mpz_t> &base, int dimension, int bit_size, char selection_char)
{
  base.resize(2 * dimension, 2 * dimension); /* Resize the matrix to a (2 * dimension) * (2 * dimension) matrix */
  /** 
   * Matrix is now filled based on specifications selected by selection_char 
   */
  if (selection_char == 'b') 
  {
  	base.gen_ntrulike2(bit_size); /* If the character is  'b', then the method first samples an integer q in space [0, 2^b-1] */ 
  }
  else 
  {
  	if (selection_char == 'q') 
  	{
  		base.gen_ntrulike2_withq(bit_size); /* If the character is 'q', then the integer is set to the provided value */	
  	}
  	else 
  	{
  		cout << "Warning: selection character MUST be 'b' or 'q'." << endl;
  	}
  	
  }
}

/**
 * Generates a q-ary matrix. It returns a 2 x 2 block matrix [[I, H], [0, q*I]],
 * where H is (d-k_param) x k_param and uniformly random modulo q. These bases correspond to
 * the SIS/LWE q-ary lattices. Goldstein-Mayer lattices correspond to k=1 and q prime
 * This functions just like 'u' method of latticegen API.
 */
void generate_qary_matrix (ZZ_mat<mpz_t> &base, int dimension, int bit_size, int k_param, char selection_char)
{
  base.resize(dimension, dimension); /* Resize the matrix to a dimension * dimension matrix */
  /** 
   * Here the matrix is created, choosing a generation method according to the selection_char input  
   */
  if (selection_char == 'b') 
  {
  	base.gen_qary(k_param, bit_size); /* If the character is  'b', then the method first samples an integer q in space [0, 2^b-1]  */
  }
  else 
  {
  	if (selection_char == 'p')
  	{
  		base.gen_qary_prime(k_param, bit_size); /* If the character is 'p', then */
  	}
  	else 
  	{
  		if (selection_char == 'q') 
  		{
  			base.gen_qary_withq(k_param, bit_size); /* If the character is 'q', then the integer is set to the provided value */	
  		}
  		else 
  		{
  			cout << "Warning: selection character MUST be 'b', 'p' or 'q'." << endl;
  		}
  	}
  }
}

/**
 * Generates a lower triangular matrix base of size dimension * dimension, where base[i][i] = 2^(dimension-i+1)^f_param for all i,
 * and base[i][j] is uniform between -base[i][j]/2 and base[i][j]/2 for all j<i. This functions exactly like 't' generation method
 * of latticegen API. 
 */
void generate_lower_triangular_matrix (ZZ_mat<mpz_t> &base, int dimension, float f_param)
{
  base.resize(dimension, dimension); /* Resize the matrix to a dimension * dimension matrix */
  base.gen_trg(f_param); /* Matrix is generated here, as described above. */
}

/**
 * Generates a lower triangular matrix base much like GenerateTriangularMatrix, except that it also receives as input a 
 * dimension-dimensional vector, read from a file. It then generates base as follows: base[i][i] =vec[i] for all i, and
 * base[i][j] is uniform betweeen -base[j][j]/2 and base[j][j]/2 for all j<i.
 */
void generate_lower_triangular_matrix_alt (ZZ_mat<mpz_t> &base, int dimension, float f_param) 
{
	FP_NR<mpfr_t> *input_vector = new FP_NR<mpfr_t>[dimension]; /* This is the input vector declaration. */
	for (int i = 0; i < dimension; i++) 
	{
		input_vector[i] = (float) i + 1;
	}
	base.resize(dimension, dimension); /* Resize the matrix to a dimension * dimension matrix */
	base.gen_trg2(input_vector); /* Matrix is generated here, as described above. */
	delete[] input_vector;
}


int main (int argc, char** argv) 
{
	int dimension = 5, bit_size = 4, bit_size2 = 4, k_param = 1, output = 0;
	char input1 = 'b', input2 = 'q', input3 = 'p';
	float f_param = 1.02;
	/* This is the lattice Base declaration. */ 
	
	/* ZZ_mat<mpz_t> integer_lattice_base: A matrix of integers. Initially, the matrix is empty. */
	
	ZZ_mat<mpz_t> integer_lattice_base;
	generate_knapsack_matrix (integer_lattice_base, dimension, bit_size);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_simdioph_matrix (integer_lattice_base, dimension, bit_size, bit_size2);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_uniformly_random_matrix (integer_lattice_base, dimension, bit_size);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_ntru_like_matrix(integer_lattice_base, dimension, bit_size, input1);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_ntru_like_matrix(integer_lattice_base, dimension, bit_size, input2);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_ntru_like_matrix_alt(integer_lattice_base, dimension, bit_size, input1);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_ntru_like_matrix_alt(integer_lattice_base, dimension, bit_size, input2);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_qary_matrix(integer_lattice_base, dimension, bit_size, k_param, input3);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_qary_matrix(integer_lattice_base, dimension, bit_size, k_param, input2);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_qary_matrix(integer_lattice_base, dimension, bit_size, k_param, input1);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_lower_triangular_matrix(integer_lattice_base, dimension, f_param);
	output = write_to_file (integer_lattice_base, "base_output");
	integer_lattice_base.clear();
	generate_lower_triangular_matrix_alt(integer_lattice_base, dimension, f_param);
	output = write_to_file (integer_lattice_base, "base_output");
	return 0;
} 
