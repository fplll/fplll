/* Copyright (C) 2015 Martin Albrecht
 *               2016 Michael Walter

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
#include <test_utils.h>

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

using namespace std;
using namespace fplll;

enum Test
{
  SVP_ENUM,
  DSVP_ENUM,
  DSVP_REDUCE
};

/**
   @brief Test if SVP function returns vector with right norm.

   @param A              input lattice
   @param b              shortest vector
   @return
*/

template <class ZT> int test_svp(ZZ_mat<ZT> &A, vector<Z_NR<mpz_t>> &b)
{
  vector<Z_NR<mpz_t>> sol_coord;   // In the LLL-reduced basis
  vector<Z_NR<mpz_t>> sol_coord2;  // In the initial basis
  vector<Z_NR<mpz_t>> solution;
  ZZ_mat<mpz_t> u;

  int status =
      lll_reduction(A, u, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);
  if (status != RED_SUCCESS)
  {
    cerr << "LLL reduction failed: " << get_red_status_str(status) << endl;
    return status;
  }

  status = shortest_vector(A, sol_coord, SVPM_PROVED, SVP_DEFAULT);

  if (status != RED_SUCCESS)
  {
    cerr << "Failure: " << get_red_status_str(status) << endl;
    return status;
  }

  vector_matrix_product(sol_coord2, sol_coord, u);
  vector_matrix_product(solution, sol_coord, A);

  Z_NR<ZT> tmp;
  Z_NR<ZT> norm_s;
  Z_NR<ZT> norm_b;

  for (int i = 0; i < A.get_cols(); i++)
  {
    tmp.mul(solution[i], solution[i]);
    norm_s.add(norm_s, tmp);

    tmp.mul(b[i], b[i]);
    norm_b.add(norm_b, tmp);
  }
  if (norm_s != norm_b)
    return 1;

  return 0;
}

/**
   @brief Compute the norm of a dual vector (specified by coefficients in the dual basis).

   @param A              input lattice
   @param b              coefficients of shortest dual vector
   @return
*/
template <class ZT>
int dual_length(FP_NR<mpfr_t> &norm, ZZ_mat<ZT> &A, const vector<Z_NR<mpz_t>> &coords)
{
  int d = coords.size();
  if (A.get_rows() != d)
  {
    cerr << "DSVP length error: Coefficient vector has wrong dimension: ";
    cerr << A.get_rows() << " vs " << d << endl;
    return 1;
  }
  vector<FP_NR<mpfr_t>> coords_d(d);
  for (int i = 0; i < d; i++)
  {
    coords_d[i] = coords[i].get_d();
  }

  ZZ_mat<mpz_t> empty_mat;
  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> gso(A, empty_mat, empty_mat, GSO_INT_GRAM);
  if (!gso.update_gso())
  {
    cerr << "GSO Failure." << endl;
    return 1;
  }
  FP_NR<mpfr_t> tmp;
  gso.get_r(tmp, d - 1, d - 1);
  tmp.pow_si(tmp, -1);

  vector<FP_NR<mpfr_t>> alpha(d);
  FP_NR<mpfr_t> mu, alpha2, r_inv;
  norm = 0.0;
  for (int i = 0; i < d; i++)
  {
    alpha[i] = coords_d[i];
    for (int j = 0; j < i; j++)
    {
      gso.get_mu(mu, i, j);
      alpha[i] -= mu * alpha[j];
    }
    gso.get_r(r_inv, i, i);
    r_inv.pow_si(r_inv, -1);
    alpha2.pow_si(alpha[i], 2);
    norm += alpha2 * r_inv;
  }

  return 0;
}

/**
   @brief Test if dual SVP function returns vector with right norm.

   @param A              input lattice
   @param b              shortest dual vector
   @return
*/

template <class ZT> int test_dual_svp(ZZ_mat<ZT> &A, vector<Z_NR<mpz_t>> &b)
{
  vector<Z_NR<mpz_t>> sol_coord;  // In the LLL-reduced basis
  vector<Z_NR<mpz_t>> solution;
  ZZ_mat<mpz_t> u;

  FP_NR<mpfr_t> normb;
  if (dual_length(normb, A, b))
  {
    return 1;
  }

  int status =
      lll_reduction(A, u, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);
  if (status != RED_SUCCESS)
  {
    cerr << "LLL reduction failed: " << get_red_status_str(status) << endl;
    return status;
  }

  status = shortest_vector(A, sol_coord, SVPM_FAST, SVP_DUAL);

  if (status != RED_SUCCESS)
  {
    cerr << "Failure: " << get_red_status_str(status) << endl;
    return status;
  }

  FP_NR<mpfr_t> norm_sol;
  if (dual_length(norm_sol, A, sol_coord))
  {
    return 1;
  }

  FP_NR<mpfr_t> error;
  error = 1;
  error.mul_2si(error, -(int)error.get_prec());
  normb += error;
  if (norm_sol > normb)
  {
    cerr << "Returned dual vector too long by more than " << error << endl;
    return 1;
  }

  return 0;
}

/**
   @brief Test if dual SVP reduction returns reduced basis.

   @param A              input lattice
   @param b              shortest dual vector
   @return
*/
template <class ZT> int test_dsvp_reduce(ZZ_mat<ZT> &A, vector<Z_NR<mpz_t>> &b)
{
  ZZ_mat<mpz_t> u;
  int d = A.get_rows();

  FP_NR<mpfr_t> normb;
  if (dual_length(normb, A, b))
  {
    return 1;
  }

  int status =
      lll_reduction(A, u, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);
  if (status != RED_SUCCESS)
  {
    cerr << "LLL reduction failed: " << get_red_status_str(status) << endl;
    return status;
  }

  ZZ_mat<mpz_t> empty_mat;
  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> gso(A, empty_mat, empty_mat, GSO_INT_GRAM);
  LLLReduction<Z_NR<mpz_t>, FP_NR<mpfr_t>> lll_obj(gso, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);

  vector<Strategy> strategies;
  BKZParam dummy(d, strategies);
  BKZReduction<Z_NR<mpz_t>, FP_NR<mpfr_t>> bkz_obj(gso, lll_obj, dummy);

  bkz_obj.svp_reduction(0, d, dummy, true);

  FP_NR<mpfr_t> norm_sol;
  Z_NR<mpz_t> zero;
  zero = 0;
  vector<Z_NR<mpz_t>> e_n(d, zero);
  e_n[d - 1] = 1;
  if (dual_length(norm_sol, A, e_n))
  {
    return 1;
  }

  FP_NR<mpfr_t> error;
  error = 1;
  error.mul_2si(error, -(int)error.get_prec());
  normb += error;
  if (norm_sol > normb)
  {
    cerr << "Last dual vector too long by more than " << error << endl;
    return 1;
  }

  return 0;
}

/**
   @brief Test if SVP function returns vector with right norm.

   @param input_filename   filename of an input lattice
   @param output_filename  filename of a shortest vector
   @return
*/

template <class ZT>
int test_filename(const char *input_filename, const char *output_filename,
                  const Test test = SVP_ENUM)
{
  ZZ_mat<ZT> A;
  int status = 0;
  status |= read_file(A, input_filename);

  vector<Z_NR<mpz_t>> b;
  status |= read_file(b, output_filename);

  switch (test)
  {
  case SVP_ENUM:
    status |= test_svp<ZT>(A, b);
    return status;
  case DSVP_ENUM:
    status |= test_dual_svp<ZT>(A, b);
    return status;
  case DSVP_REDUCE:
    status |= test_dsvp_reduce<ZT>(A, b);
    return status;
  }

  cerr << "Unknown test." << endl;
  return 1;
}

/**
   @brief Run SVP tests.

   @return
*/

int main()
{

  int status = 0;
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_svp_in",
                                 TESTDATADIR "/tests/lattices/example_svp_out");
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_dsvp_in",
                                 TESTDATADIR "/tests/lattices/example_dsvp_out", DSVP_ENUM);
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_dsvp_in",
                                 TESTDATADIR "/tests/lattices/example_dsvp_out", DSVP_REDUCE);

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
