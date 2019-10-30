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
