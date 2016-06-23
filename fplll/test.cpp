/* Copyright (C) 2011 Xavier Pujol.

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


#include <fplll.h>

using namespace std;
using namespace fplll;

int main(void) {
  mpz_t q;
  ZZ_mat<mpz_t> M(3, 3);
  ZZ_mat<mpz_t> U(3, 3);
  M.gen_uniform(3);
  mpz_init_set_str(q, "FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831", 16); 
  M[2][2].set(q);
  cout << M << endl;
  lllReduction(M, U, 0.99, 0.51, LM_WRAPPER);
  cout << M << endl;
  cout << U << endl;
  mpz_clear(q);
  return 0;
}
