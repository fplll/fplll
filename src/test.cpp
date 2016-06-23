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
#define N 56


  Pruner<FP_NR<double>> pru;

  double rs[N];
  double dcost[N];
  double pr[N];
  


  for (int i = 0; i < N; ++i)
  {
    rs[i] = pow(1.06, - i);
    pr[i] = 1.;
  }


  pru.enumeration_radius = .85;
  pru.target_success_proba = .50;
  pru.preproc_cost = 1e10;

  pru.load_basis_shape(N, rs);

  cerr << "un-pruned cost" << pru.get_cost(pr) << endl;


  pru.optimize_pruning_coeffs(pr);
  cerr << "Success Proba " << pru.get_svp_success_proba(pr) << endl;
  cerr << "Cost " << pru.get_cost(pr) << endl;


  for (int i = 0; i < N; ++i)
  {
    cerr << pr[i] << ", ";
  }
  cerr << endl;


  return 0;
}

