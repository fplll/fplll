/* Copyright (C) 2016 Martin Albrecht

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

#include "io/json.hpp"
#include <cstring>
#include <fplll.h>
#include <test_utils.h>

using json = nlohmann::json;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

using namespace std;
using namespace fplll;

/**
   @brief Test BKZ reduction.

   @param A                test matrix
   @param block_size       block size
   @param float_type       floating point type to test
   @param flags            flags to use
   @param prec             precision if mpfr is used

   @return zero on success.
*/

template <class ZT>
int test_bkz(ZZ_mat<ZT> &A, const int block_size, FloatType float_type, int flags = BKZ_DEFAULT,
             int prec = 0)
{

  int status = 0;

  // zero on success
  status = bkz_reduction(A, block_size, flags, float_type, prec);
  if (status != RED_SUCCESS)
  {
    cerr << "BKZ reduction failed with error '" << get_red_status_str(status);
    cerr << " for float type " << FLOAT_TYPE_STR[float_type] << endl;
  }
  return status;
}

/**
   @brief Test BKZ strategy interface.

   @param A                test matrix
   @param block_size       block size

   @return zero on success.
*/

template <class ZT>
int test_bkz_param(ZZ_mat<ZT> &A, const int block_size, int flags = BKZ_DEFAULT,
                   string dump_gso_filename = string())
{

  int status = 0;

  vector<Strategy> strategies;
  for (long b = 0; b <= block_size; b++)
  {
    Strategy strategy = Strategy::EmptyStrategy(b);
    if (b == 10)
    {
      strategy.preprocessing_block_sizes.emplace_back(5);
    }
    else if (b == 20)
    {
      strategy.preprocessing_block_sizes.emplace_back(10);
    }
    else if (b == 30)
    {
      strategy.preprocessing_block_sizes.emplace_back(15);
    }
    strategies.emplace_back(std::move(strategy));
  }

  BKZParam params(block_size, strategies);
  params.flags             = flags;
  params.dump_gso_filename = dump_gso_filename;
  // zero on success
  status = bkz_reduction(&A, NULL, params, FT_DEFAULT, 53);
  if (status != RED_SUCCESS)
  {
    cerr << "BKZ parameter test failed with error '" << get_red_status_str(status) << "'" << endl;
  }
  return status;
}

/**
   @brief Test BKZ with pruning.

   @param A                test matrix
   @param block_size       block size

   @return zero on success.
*/

template <class ZT>
int test_bkz_param_linear_pruning(ZZ_mat<ZT> &A, const int block_size, int flags = BKZ_DEFAULT)
{

  int status = 0;
  vector<Strategy> strategies;
  for (long b = 0; b < block_size; b++)
  {
    Strategy strategy = Strategy::EmptyStrategy(b);
    if (b == 10)
    {
      strategy.preprocessing_block_sizes.emplace_back(5);
    }
    else if (b == 20)
    {
      strategy.preprocessing_block_sizes.emplace_back(10);
    }
    else if (b == 30)
    {
      strategy.preprocessing_block_sizes.emplace_back(15);
    }
    strategies.emplace_back(std::move(strategy));
  }

  Strategy strategy;
  strategy.pruning_parameters.emplace_back(
      PruningParams::LinearPruningParams(block_size, block_size / 2));
  strategies.emplace_back(std::move(strategy));

  BKZParam params(block_size, strategies);
  params.flags = flags;
  // zero on success
  status = bkz_reduction(&A, NULL, params, FT_DEFAULT, 53);
  if (status != RED_SUCCESS)
  {
    cerr << "BKZ parameter test failed with error '" << get_red_status_str(status) << "'" << endl;
  }
  return status;
}

template <class ZT>
int test_bkz_param_pruning(ZZ_mat<ZT> &A, const int block_size, int flags = BKZ_DEFAULT)
{

  int status = 0;

  vector<Strategy> strategies = load_strategies_json(TESTDATADIR "/strategies/default.json");
  BKZParam params(block_size, strategies);
  params.flags = flags;
  // zero on success
  status = bkz_reduction(&A, NULL, params, FT_DEFAULT, 53);
  if (status != RED_SUCCESS)
  {
    cerr << "BKZ parameter test failed with error '" << get_red_status_str(status) << "'" << endl;
  }
  return status;
}

/**
   @brief Test BKZ_DUMP_GSO for a matrix d × (d+1) integer relations matrix with bit size b (copied
   from test_int_rel)

   @param d                dimension
   @param b                bit size
   @param block_size       block size
   @param float_type       floating point type to test
   @param flags            flags to use

   @return zero on success.
 */
template <class ZT>
int test_int_rel_bkz_dump_gso(int d, int b, const int block_size,
                              int flags = BKZ_DEFAULT | BKZ_DUMP_GSO)
{
  ZZ_mat<ZT> A, B;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  B          = A;
  int status = 0;
  // TODO: maybe not safe.
  string file_bkz_dump_gso = tmpnam(nullptr);
  status |= test_bkz_param<ZT>(B, block_size, flags, file_bkz_dump_gso);

  if (status != 0)
  {
    cerr << "Error in test_bkz_param." << endl;
    return status;
  }

  json js;
  std::ifstream fs(file_bkz_dump_gso);
  if (fs.fail())
  {
    cerr << "File " << file_bkz_dump_gso << " cannot be loaded." << endl;
    return 1;
  }
  fs >> js;
  int loop    = -1;
  double time = 0.0;
  for (auto i : js)
  {
    // Verify if there are as much norms as there are rows in A
    if (A.get_rows() != (int)i["norms"].size())
    {
      cerr << "The array \"norms\" does not contain enough values (" << A.get_rows()
           << " expected but " << i["norms"].size() << " found)." << endl;
      return 1;
    }

    // Extract data from json file
    const string step_js = i["step"];
    const int loop_js    = i["loop"];
    const double time_js = i["time"];
    // Verify if loop of Input and Output have loop = -1
    if (step_js.compare("Input") == 0 || step_js.compare("Output") == 0)
    {
      if (loop_js != -1)
      {
        cerr << "Steps Input or Output are not with \"loop\" = -1." << endl;
        return 1;
      }
    }
    else
    {
      // Verify that loop increases
      loop++;
      if (loop_js != loop)
      {
        cerr << "Loop does not increase." << endl;
        return 1;
      }
      // Verify that time increases
      if (time > time_js)
      {
        cerr << "Time does not increase." << endl;
        return 1;
      }
      time = time_js;
    }
  }

  return 0;
}

/**
   @brief Test BKZ for matrix stored in file pointed to by `input_filename`.

   @param input_filename   a path
   @param block_size       block size
   @param float_type       floating point type to test
   @param flags            flags to use
   @param prec             precision if mpfr is used

   @return zero on success
*/

template <class ZT>
int test_filename(const char *input_filename, const int block_size,
                  FloatType float_type = FT_DEFAULT, int flags = BKZ_DEFAULT, int prec = 0)
{
  ZZ_mat<ZT> A, B;
  int status = 0;
  status |= read_matrix(A, input_filename);
  B = A;
  status |= test_bkz<ZT>(A, block_size, float_type, flags, prec);
  status |= test_bkz_param<ZT>(B, block_size);
  return status;
}

/**
   @brief Construct d × (d+1) integer relations matrix with bit size b and test BKZ.

   @param d                dimension
   @param b                bit size
   @param block_size       block size
   @param float_type       floating point type to test
   @param flags            flags to use
   @param prec             precision if mpfr is used

   @return zero on success
*/

template <class ZT>
int test_int_rel(int d, int b, const int block_size, FloatType float_type = FT_DEFAULT,
                 int flags = BKZ_DEFAULT, int prec = 0)
{
  ZZ_mat<ZT> A, B;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  B          = A;
  int status = 0;
  status |= test_bkz<ZT>(A, block_size, float_type, flags | BKZ_VERBOSE, prec);
  status |= test_bkz_param<ZT>(B, block_size);
  status |= test_bkz_param_linear_pruning<ZT>(B, block_size);
  status |= test_bkz_param_pruning<ZT>(B, block_size);
  return status;
}

int test_linear_dep()
{
  ZZ_mat<mpz_t> A;
  std::stringstream("[[1 2 3]\n [4 5 6]\n [7 8 9]]\n") >> A;
  return test_bkz_param<mpz_t>(A, 3);
}

int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;

  status |= test_linear_dep();
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/dim55_in", 10, FT_DEFAULT,
                                 BKZ_DEFAULT | BKZ_AUTO_ABORT);
#ifdef FPLLL_WITH_QD
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/dim55_in", 10, FT_DD,
                                 BKZ_SD_VARIANT | BKZ_AUTO_ABORT);
#endif
  status |=
      test_filename<mpz_t>(TESTDATADIR "/tests/lattices/dim55_in", 10, FT_DEFAULT, BKZ_SLD_RED);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/dim55_in", 20, FT_MPFR,
                                 BKZ_DEFAULT | BKZ_AUTO_ABORT, 128);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/dim55_in", 20, FT_MPFR,
                                 BKZ_SD_VARIANT | BKZ_AUTO_ABORT, 128);
  status |=
      test_filename<mpz_t>(TESTDATADIR "/tests/lattices/dim55_in", 20, FT_MPFR, BKZ_SLD_RED, 128);

  status |= test_int_rel<mpz_t>(50, 1000, 10, FT_DOUBLE, BKZ_DEFAULT | BKZ_AUTO_ABORT);
  status |= test_int_rel<mpz_t>(50, 1000, 10, FT_DOUBLE, BKZ_SD_VARIANT | BKZ_AUTO_ABORT);
  status |= test_int_rel<mpz_t>(50, 1000, 10, FT_DOUBLE, BKZ_SLD_RED);
  status |= test_int_rel<mpz_t>(50, 1000, 15, FT_MPFR, BKZ_DEFAULT | BKZ_AUTO_ABORT, 100);
  status |= test_int_rel<mpz_t>(50, 1000, 15, FT_MPFR, BKZ_SD_VARIANT | BKZ_AUTO_ABORT, 100);
  status |= test_int_rel<mpz_t>(50, 1000, 15, FT_MPFR, BKZ_SLD_RED, 100);

  status |= test_int_rel<mpz_t>(30, 2000, 10, FT_DPE, BKZ_DEFAULT | BKZ_AUTO_ABORT);
  status |= test_int_rel<mpz_t>(30, 2000, 10, FT_DPE, BKZ_SD_VARIANT | BKZ_AUTO_ABORT);
  status |= test_int_rel<mpz_t>(30, 2000, 10, FT_DPE, BKZ_SLD_RED);
  status |= test_int_rel<mpz_t>(30, 2000, 10, FT_MPFR, BKZ_DEFAULT | BKZ_AUTO_ABORT, 53);
  status |= test_int_rel<mpz_t>(30, 2000, 10, FT_MPFR, BKZ_SD_VARIANT | BKZ_AUTO_ABORT, 53);
  status |= test_int_rel<mpz_t>(30, 2000, 10, FT_MPFR, BKZ_SLD_RED, 53);

  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_DEFAULT,
                                 BKZ_SD_VARIANT);
  status |=
      test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_DEFAULT, BKZ_SLD_RED);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_DOUBLE);
  status |=
      test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_DOUBLE, BKZ_SD_VARIANT);
  status |=
      test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_DOUBLE, BKZ_SLD_RED);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_MPFR,
                                 BKZ_AUTO_ABORT, 212);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_MPFR,
                                 BKZ_SD_VARIANT | BKZ_AUTO_ABORT, 212);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_MPFR,
                                 BKZ_SLD_RED | BKZ_AUTO_ABORT, 212);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_DOUBLE);
  status |=
      test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_DOUBLE, BKZ_SD_VARIANT);
  status |=
      test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_in", 10, FT_DOUBLE, BKZ_SLD_RED);

  // Test BKZ_DUMP_GSO
  status |= test_int_rel_bkz_dump_gso<mpz_t>(50, 1000, 15, BKZ_DEFAULT | BKZ_DUMP_GSO);

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
