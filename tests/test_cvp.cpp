/* Copyright (C) 2016 Guillaume Bonnoron

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

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif
using namespace std;
using namespace fplll;

/**
   @brief Read matrix from `input_filename`.

   @param A
   @param input_filename
   @return zero if A is not empty.
*/

template <class ZT> int read_matrix(ZZ_mat<ZT> &A, const char *input_filename)
{
  istream *is = new ifstream(input_filename);
  *is >> A;
  int status = 0;
  if (A.empty())
  {
    status = 1;
    cerr << "File " << input_filename << " was (probably) not opened." << endl;
  }
  delete is;
  return status;
}

/**
   @brief Read vector from `input_filename` into `b`.

   @param b                vector
   @param input_filename   filename
   @return zero if b is not empty.
*/

template <class ZT> int read_vector(vector<Z_NR<ZT>> &b, const char *input_filename)
{
  istream *is = new ifstream(input_filename);
  *is >> b;
  int status = 0;
  if (b.empty())
  {
    status = 1;
    cerr << "File " << input_filename << " was (probably) not opened." << endl;
  }
  delete is;
  return status;
}

/**
   @brief Test if CVP function returns correct vector.

   @param A              input lattice
   @param b              expected vector
   @return
*/

template <class ZT> int test_cvp(ZZ_mat<ZT> &A, IntVect &target, IntVect &b, const int method)
{
  IntVect sol_coord;  // In the LLL-reduced basis
  IntVect solution;
  IntMatrix u;

  // cerr << "A " << endl << A << endl;
  int status =
      lll_reduction(A, u, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);
  if (status != RED_SUCCESS)
  {
    cerr << "LLL reduction failed: " << get_red_status_str(status) << endl;
    return status;
  }

  status = closest_vector(A, target, sol_coord, method);

  if (status != RED_SUCCESS)
  {
    cerr << "Failure: " << get_red_status_str(status) << endl;
    return status;
  }

  vector_matrix_product(solution, sol_coord, A);

  // cerr << "A red" << endl << A << endl;
  // cerr << "sol_coord : " << sol_coord << endl;
  // cerr << "solution : " << solution << endl;
  // cerr << "expected : " << b << endl;
  // cerr << "target : " << target << endl;

  bool correct = true;
  for (int i = 0; i < A.get_cols(); i++)
  {
    correct = correct && (solution[i] == b[i]);
  }
  if (!correct)
    return 1;

  return 0;
}

/**
   @brief Test if CVP function returns the correct vector

   @param input_filename_lattice   filename of an input lattice
   @param input_filename_target    filename of a target vector
   @param output_filename  filename of the expected vector
   @return
*/

template <class ZT>
int test_filename(const char *input_filename_lattice, const char *input_filename_target,
                  const char *output_filename, const int method = CVPM_FAST)
{
  ZZ_mat<ZT> A;
  int status = 0;
  status |= read_matrix(A, input_filename_lattice);

  IntVect t;
  status |= read_vector(t, input_filename_target);

  IntVect b;
  status |= read_vector(b, output_filename);

  status |= test_cvp<ZT>(A, t, b, method);
  return status;
}

/**
   @brief Run CVP tests.

   @param argc             ignored
   @param argv             ignored
   @return
*/

int main(int argc, char *argv[])
{

  int status = 0;
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice",
                                 TESTDATADIR "/tests/lattices/example_cvp_in_target",
                                 TESTDATADIR "/tests/lattices/example_cvp_out");
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2",
                                 TESTDATADIR "/tests/lattices/example_cvp_in_target2",
                                 TESTDATADIR "/tests/lattices/example_cvp_out2");
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3",
                                 TESTDATADIR "/tests/lattices/example_cvp_in_target3",
                                 TESTDATADIR "/tests/lattices/example_cvp_out3");
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4",
                                 TESTDATADIR "/tests/lattices/example_cvp_in_target4",
                                 TESTDATADIR "/tests/lattices/example_cvp_out4", CVPM_PROVED);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5",
                                 TESTDATADIR "/tests/lattices/example_cvp_in_target5",
                                 TESTDATADIR "/tests/lattices/example_cvp_out5", CVPM_PROVED);

  if (status == 0)
  {
    cerr << "All tests passed." << endl;
    return 0;
  }
  else
  {
    return -1;
  }

  return 0;
}
