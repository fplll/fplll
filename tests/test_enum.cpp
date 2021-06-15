/* Copyright (C) 2019 Martin Albrecht

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
#include <fplll/fplll.h>
#include <test_utils.h>

using namespace fplll;

template <class FT> int test_enum(size_t d)
{
  RandGen::init_with_seed(0x1337);
  ZZ_mat<mpz_t> A = ZZ_mat<mpz_t>(100, 100);
  A.gen_qary_withq(50, 7681);
  lll_reduction(A);
  ZZ_mat<mpz_t> U;
  MatGSO<Z_NR<mpz_t>, FP_NR<FT>> M(A, U, U, 0);
  M.update_gso();

  FastEvaluator<FP_NR<FT>> evaluator;
  Enumeration<Z_NR<mpz_t>, FP_NR<FT>> enum_obj(M, evaluator);
  FP_NR<FT> max_dist;
  M.get_r(max_dist, 0, 0);
  max_dist *= 0.99;
  enum_obj.enumerate(0, d, max_dist, 0);
  if (evaluator.empty())
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

/**
   @brief Test if list_CVP via enumeration function returns the correct amount of vectors
   @return
*/

template <class FT> int test_list_cvp()
{
  ZZ_mat<mpz_t> u;
  int status = 0;
  ZZ_mat<mpz_t> A;
  status |= read_file(A, TESTDATADIR "tests/lattices/example_list_cvp_in_lattice");

  status |= lll_reduction(A);
  if (status != RED_SUCCESS)
  {
    cerr << "LLL reduction failed: " << get_red_status_str(status) << endl;
    return status;
  }

  // Search for up to 999999 vectors, up to radius 32.5 around the origin
  // the right answer is 196561
  FP_NR<FT> rad                   = 32.5;
  FP_NR<FT> half_rad              = 16.5;
  const unsigned int right_answer = 196561;

  // Tests with two targets: 0, and something very close to 0.

  // HOLE: Not sure how to set that up
  size_t d = A.get_rows();
  if (d != 24)
  {
    cerr << "Expected a lattice of dimension 24, got : " << d << endl;
    return 1;
  }

  ZZ_mat<mpz_t> empty_mat;
  MatGSO<Z_NR<mpz_t>, FP_NR<FT>> gso(A, empty_mat, empty_mat, GSO_INT_GRAM);

  {
    gso.update_gso();
    FastEvaluator<FP_NR<FT>> evaluator(999999);
    Enumeration<Z_NR<mpz_t>, FP_NR<FT>> enum_obj(gso, evaluator);

    std::vector<FP_NR<FT>> target(d, 0.0);
    enum_obj.enumerate(0, d, rad, 0, target);
    if (evaluator.size() != right_answer)
    {
      cerr << "list CVP failed (center at 0), expected 196561 solutions, got : " << evaluator.size()
           << endl;
      return 1;
    }
    // cerr << "list CVP (at 0) PASSED : " << evaluator.size() << endl;
  }

  {
    gso.update_gso();
    FastEvaluator<FP_NR<FT>> evaluator(999999);
    Enumeration<Z_NR<mpz_t>, FP_NR<FT>> enum_obj(gso, evaluator);

    std::vector<FP_NR<FT>> target(d, 0.0001);
    enum_obj.enumerate(0, d, rad, 0, target);
    if (evaluator.size() != right_answer)
    {
      cerr << "list CVP failed (center near 0), expected 196561 solutions, got : "
           << evaluator.size() << endl;
      return 1;
    }
    // cerr << "list CVP (near 0) PASSED : " << evaluator.size() << endl;
  }

  for (size_t rep = 0; rep < 24; ++rep)
  {
    FastEvaluator<FP_NR<FT>> evaluator(999999);
    Enumeration<Z_NR<mpz_t>, FP_NR<FT>> enum_obj(gso, evaluator);

    std::vector<FP_NR<FT>> can_target(d, 0.0);

    // Generate a point in the half lattice
    for (size_t i = 0; i < d; ++i)
    {
      size_t c = (i == rep || rand() % 2);  // make sure at least one of them is non-zero
      // std::cerr << c;
      if (!c)
        continue;
      for (size_t j = 0; j < d; ++j)
      {
        can_target[j] += 0.5 * A[i][j].get_d();
      }
    }

    std::vector<FP_NR<FT>> target(d, 0.0);
    // Convert it to GSO basis
    for (size_t i = 0; i < d; ++i)
    {
      for (size_t j = 0; j < d; ++j)
      {
        target[i] += A[i][j].get_d() * can_target[j];
      }

      for (size_t j = 0; j < i; ++j)
      {
        FP_NR<FT> mu_ij;
        gso.get_mu(mu_ij, i, j);
        target[i] -= target[j] * mu_ij;
      }
    }
    for (size_t i = 0; i < d; ++i)
    {
      FP_NR<FT> r_ii;
      gso.get_r(r_ii, i, i);
      target[i] /= r_ii;
    }

    enum_obj.enumerate(0, d, half_rad, 0, target);

    if (evaluator.size() != 2 && evaluator.size() != 48)
    {
      cerr << "list CVP failed (halfway), expected 2 or 48 solutions, got : " << evaluator.size()
           << endl;
      return 1;
    }
    // cerr << "list CVP failed (halfway) PASSED : " << evaluator.size() << endl;
  }

  return 0;
}

bool callback_firstf(size_t n, enumf *new_sol_coord, void *ctx)
{
  if (new_sol_coord[0] == static_cast<double *>(ctx)[0])
  {
    return true;
  }
  return false;
}

template <class FT> int test_callback_enum(size_t d)
{
  RandGen::init_with_seed(0x1337);
  ZZ_mat<mpz_t> A = ZZ_mat<mpz_t>(100, 100);
  A.gen_qary_withq(50, 7681);
  lll_reduction(A);
  ZZ_mat<mpz_t> U;
  MatGSO<Z_NR<mpz_t>, FP_NR<FT>> M(A, U, U, 0);
  M.update_gso();

  enumf ctx = 2;
  CallbackEvaluator<FP_NR<FT>> evaluator(callback_firstf, &ctx);
  Enumeration<Z_NR<mpz_t>, FP_NR<FT>> enum_obj(M, evaluator);
  FP_NR<FT> max_dist;
  M.get_r(max_dist, 0, 0);
  max_dist *= 0.99;
  enum_obj.enumerate(0, d, max_dist, 0);
  if (evaluator.empty())
  {
    return 1;
  }
  else
  {
    if (evaluator.begin()->second[0].get_si() == 2)
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }
}

int main(int argc, char *argv[])
{
  int status = 0;
  status |= test_enum<double>(30);
  status |= test_callback_enum<double>(40);
  status |= test_list_cvp<double>();
  if (status == 0)
  {
    std::cerr << "All tests passed." << std::endl;
    return 0;
  }
  else
  {
    return -1;
  }
}
