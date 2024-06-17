/* Copyright (C) 2016 Martin Albrecht
   Copyright (C) 2019 Koen de Boer & Wessel van Woerden

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

#include "test_utils.h"
#include <cstring>
#include <fplll/fplll.h>
#include <fplll/gso_gram.h>
#include <fplll/gso_interface.h>
#include <fplll/io/json.hpp>

using json = nlohmann::json;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

using namespace std;
using namespace fplll;

template <class ZT> bool gram_is_equal(ZZ_mat<ZT> b, ZZ_mat<ZT> G)
{
  int r = b.r;
  int c = b.c;
  ZZ_mat<ZT> G_reduced;
  G_reduced.resize(r, r);
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < r; j++)
    {
      (b)[i].dot_product(G_reduced(i, j), (b)[j], c);
    }
  }
  // ------------------------------------------------
  // ************************************************

  // _______________________________________________
  // -----------------------------------------------
  // Test whether G_reduced = G
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < i; j++)
    {
      if (G(i, j) != G_reduced(i, j))
      {
        cerr << "The gram-representation and the basis-representation of the same lattice have an "
                "unequal gram matrix.\n";
        return 1;
      }
    }
  }
  return 0;
}

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
  FloatType sel_ft = (float_type != FT_DEFAULT) ? float_type : FT_DOUBLE;

  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  /* lllwrapper (no FloatType needed, -m ignored) */
  if (flags & BKZ_NO_LLL)
    zeros_last(A, U, UT);
  else
  {
    Wrapper wrapper(A, U, UT, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
    if (!wrapper.lll())
      return wrapper.status;
  }

  // _______________________________________________
  // -----------------------------------------------
  // Create the Gram matrix G of the basis A

  ZZ_mat<ZT> G;
  int r = A.r;
  int c = A.c;
  G.resize(r, r);
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < r; j++)
    {
      A[i].dot_product(G(i, j), A[j], c);
    }
  }
  // ------------------------------------------------
  // ************************************************

  // _______________________________________________
  // -----------------------------------------------
  // Create a MatGSO-object M (basis gso) for A
  // and a MatGSOGram-object Mgram (gram gso) for G.

  vector<Strategy> strategies;
  BKZParam param(block_size, strategies);
  param.flags = flags;

  if (sel_ft == FT_DOUBLE)
  {
    MatGSO<Z_NR<ZT>, FP_NR<double>> M(A, U, UT, 0);
    M.update_gso();
    MatGSOGram<Z_NR<ZT>, FP_NR<double>> Mgram(G, U, UT, 1);
    Mgram.update_gso();

    LLLReduction<Z_NR<ZT>, FP_NR<double>> LLLObj(M, LLL_DEF_DELTA, LLL_DEF_ETA, 0);
    BKZReduction<Z_NR<ZT>, FP_NR<double>> BKZObj(M, LLLObj, param);
    BKZObj.bkz();

    LLLReduction<Z_NR<ZT>, FP_NR<double>> LLLObjgram(Mgram, LLL_DEF_DELTA, LLL_DEF_ETA, 0);
    BKZReduction<Z_NR<ZT>, FP_NR<double>> BKZObjgram(Mgram, LLLObjgram, param);
    BKZObjgram.bkz();

    return gram_is_equal(A, G);
  }
  else if (sel_ft == FT_MPFR)
  {
    int old_prec = FP_NR<mpfr_t>::set_prec(prec);

    MatGSO<Z_NR<ZT>, FP_NR<mpfr_t>> M(A, U, UT, 1);
    M.update_gso();
    MatGSOGram<Z_NR<ZT>, FP_NR<mpfr_t>> Mgram(G, U, UT, 1);
    Mgram.update_gso();

    LLLReduction<Z_NR<ZT>, FP_NR<mpfr_t>> LLLObj(M, LLL_DEF_DELTA, LLL_DEF_ETA, 0);
    BKZReduction<Z_NR<ZT>, FP_NR<mpfr_t>> BKZObj(M, LLLObj, param);
    BKZObj.bkz();

    LLLReduction<Z_NR<ZT>, FP_NR<mpfr_t>> LLLObjgram(Mgram, LLL_DEF_DELTA, LLL_DEF_ETA, 0);
    BKZReduction<Z_NR<ZT>, FP_NR<mpfr_t>> BKZObjgram(Mgram, LLLObjgram, param);
    BKZObjgram.bkz();

    FP_NR<mpfr_t>::set_prec(old_prec);

    return gram_is_equal(A, G);
  }
  else
  {
    cerr << "Type not supported for test" << endl;
    return 0;
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
  printf("%s %d\n", input_filename, block_size);
  ZZ_mat<ZT> A, B;
  int status = 0;
  status |= read_file(A, input_filename);
  B = A;
  status |= test_bkz<ZT>(A, block_size, float_type, flags, prec);
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
  cerr << "test_int_rel " << d << "  " << b << " " << block_size << endl;
  ZZ_mat<ZT> A, B;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  B          = A;
  int status = 0;
  status |= test_bkz<ZT>(A, block_size, float_type, flags | BKZ_VERBOSE, prec);
  return status;
}

int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;

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

  // status |= test_int_rel<mpz_t>(30, 2000, 10, FT_DPE, BKZ_DEFAULT | BKZ_AUTO_ABORT);
  // status |= test_int_rel<mpz_t>(30, 2000, 10, FT_DPE, BKZ_SD_VARIANT | BKZ_AUTO_ABORT);
  // status |= test_int_rel<mpz_t>(30, 2000, 10, FT_DPE, BKZ_SLD_RED);
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
