/* Copyright (C) 2015 Martin Albrecht

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
#include <gso.h>
#include <gso_gram.h>
#include <gso_interface.h>
#include <lll.h>
#include <test_utils.h>

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

template <class ZT, class FT> int is_already_reduced(ZZ_mat<ZT> &A, Matrix<Z_NR<ZT>> &G)
{
  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  MatGSO<Z_NR<ZT>, FP_NR<FT>> M(A, U, UT, 1);
  MatGSOGram<Z_NR<ZT>, FP_NR<FT>> M2(G, U, UT, 1);

  int is_reduced  = is_lll_reduced<Z_NR<ZT>, FP_NR<FT>>(M, LLL_DEF_DELTA, LLL_DEF_ETA);
  int is_greduced = is_lll_reduced<Z_NR<ZT>, FP_NR<FT>>(M2, LLL_DEF_DELTA, LLL_DEF_ETA);
  if (is_reduced)
    cerr << "is_lll_reduced reports success when it should not" << endl;
  if (is_greduced)
    cerr << "is_lll_greduced reports success when it should not" << endl;
  return (is_reduced || is_greduced);
}

/**
   @brief Tests whether LLL applied on A is equivalent to LLLGram applied on G = A*A^T

   @param A basis matrix of a lattice
   @return zero on success
*/

template <class ZT, class FT> int test_lll(ZZ_mat<ZT> &A)
{

  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

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

  MatGSO<Z_NR<ZT>, FP_NR<FT>> M(A, U, UT, 1);
  M.update_gso();
  MatGSOGram<Z_NR<ZT>, FP_NR<FT>> Mgram(G, U, UT, 1);
  Mgram.update_gso();
  // ------------------------------------------------
  // ************************************************

  // _________________________________________________
  // -------------------------------------------------
  // Test whether A and G are not already LLL-reduced.
  if (is_already_reduced<ZT, FT>(A, G))
  {
    cerr << "The matrices are already LLL-reduced";
    return 1;
  }
  // ------------------------------------------------
  // ************************************************

  // _________________________________________________
  // -------------------------------------------------
  // Make two LLLObjects. One of the basis MatGSO-object M
  // one of the gram MatGSOGram-object Mgram.
  LLLReduction<Z_NR<ZT>, FP_NR<FT>> LLLObj(M, LLL_DEF_DELTA, LLL_DEF_ETA, 0);
  LLLReduction<Z_NR<ZT>, FP_NR<FT>> LLLObjgram(Mgram, LLL_DEF_DELTA, LLL_DEF_ETA, 0);

  // and LLL reduce both objects
  LLLObj.lll();
  LLLObjgram.lll();
  // ------------------------------------------------
  // ************************************************

  // _________________________________________________
  // -------------------------------------------------
  // Check whether M and Mgram are really reduced after LLL reduction
  int is_reduced  = is_lll_reduced<Z_NR<ZT>, FP_NR<FT>>(M, LLL_DEF_DELTA, LLL_DEF_ETA);
  int is_greduced = is_lll_reduced<Z_NR<ZT>, FP_NR<FT>>(Mgram, LLL_DEF_DELTA, LLL_DEF_ETA);

  if (is_reduced != 1 || is_greduced != 1)
  {
    if (is_reduced != 1)
    {
      cerr << "The basis GSO-object is not LLL-reduced after calling LLL\n";
    }
    if (is_greduced != 1)
    {
      cerr << "The gram GSO-object is not LLL-reduced after calling LLL\n";
    }
    return 1;
  }
  // ------------------------------------------------
  // ************************************************

  //
  // After this reduction, A (the input matrix) and G (gram matrix)
  // are LLL-reduced accordingly. Therefore we should check whether
  // A *A^T = G
  //

  // _______________________________________________
  // -----------------------------------------------
  // Create the Gram matrix G_reduced of the basis A_reduced

  ZZ_mat<ZT> G_reduced;
  G_reduced.resize(r, r);
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < r; j++)
    {
      (M.b)[i].dot_product(G_reduced(i, j), (M.b)[j], c);
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

  // Return 0 on success.
  return 0;
}

template <class ZT, class FT> int test_filename(const char *input_filename)
{
  ZZ_mat<ZT> A;
  int status = 0;
  status |= read_matrix(A, input_filename);
  status |= test_lll<ZT, FT>(A);
  return status;
}

/**
   @brief Construct d × (d+1) integer relations matrix with bit size b and test LLL.

   @param d                dimension
   @param b                bit size

   @return zero on success
*/

template <class ZT, class FT> int test_int_rel(int d, int b)
{
  ZZ_mat<ZT> A;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  return test_lll<ZT, FT>(A);
}

int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;

  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, double>(50, 20);
  status |= test_int_rel<mpz_t, double>(40, 10);

  status |= test_int_rel<mpz_t, double>(100, 200);

  status |= test_int_rel<mpz_t, double>(200, 100);  

  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, mpfr_t>(50, 20);
  status |= test_int_rel<mpz_t, mpfr_t>(40, 10);

#ifdef FPLLL_WITH_LONG_DOUBLE
  status |= test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |=
      test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |=
      test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |=
      test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |=
      test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, long double>(50, 20);
  status |= test_int_rel<mpz_t, long double>(40, 10);

#endif
#ifdef FPLLL_WITH_QD
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, dd_real>(50, 20);
  status |= test_int_rel<mpz_t, dd_real>(40, 10);

  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, qd_real>(50, 20);
  status |= test_int_rel<mpz_t, qd_real>(40, 10);

#endif
#ifdef FPLLL_WITH_DPE
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, dpe_t>(50, 20);
  status |= test_int_rel<mpz_t, dpe_t>(40, 10);

#endif

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
