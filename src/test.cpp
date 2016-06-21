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
  ZZ_mat<mpz_t> M(10, 10);
  ZZ_mat<mpz_t> U(0, 0);
  M.gen_uniform(3);
  //mpz_init_set_str(q, "831", 16); 
  M[2][2].set(q);
  cout << M << endl;
  MatGSO<Z_NR<mpz_t>, FP_NR<double> > gso(M,U,U,GSO_DEFAULT);
  gso.updateGSO();


  cout << M << endl;
  mpz_clear(q);
  Pruner<FP_NR<double>> pru;
  pru.load_basis_shape<Z_NR<mpz_t>,FP_NR<double>>(gso);
  double pr[10];

  for (int i = 0; i < 10; ++i)
  {
    pr[i] = 1.;
  }

  cerr << pru.get_svp_success_proba(pr) << endl;

  for (int i = 0; i < 10; ++i)
  {
    pr[i] = .5;
  }

  cerr << pru.get_svp_success_proba(pr) << endl;


  return 0;
}
