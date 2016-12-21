/* Copyright (C) 2016 Leo Ducas

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

#ifdef FPLLL_WITH_QD
#include <qd/dd_real.h>
#endif

using namespace std;
using namespace fplll;

/**
   @brief Read matrix from `input_filename`.

   @param A
   @param input_filename
   @return
*/

#define N 56
#define D N / 2
#define Nbis 24
#define Dbis 12

template <class FT> class Pruner<FT>::TestPruner
{
public:
  Pruner<FT> pru;
  int test_enforce()
  {
    vector<double> rs;
    int status = 0;
    for (int i = 0; i < N; ++i)
    {
      rs.emplace_back(1.);
    }
    pru.load_basis_shape(rs);
    Pruner<FT>::evec b;
    for (int i = 0; i < D; ++i)
    {
      b[i] = .3 + (1. * ((i * i * i) % D)) / (2 * D);
    }
    int j = D / 2;
    b[j]  = .5;

    FT old_bj = b[j];
    pru.enforce_bounds(b, j);

    status |= !((b[j] >= old_bj) && (b[j] <= old_bj));
    status |= !((b[D - 1] >= 1.) && (b[D - 1] <= 1.));
    cerr << "Checking enforced bounds: Increasing, 1-terminated, keeping b[13] = .5" << endl;
    for (int i = 0; i < D - 1; ++i)
    {
      status |= !(b[i + 1] >= b[i]);
      cerr << b[i] << " , ";
    }
    cerr << b[D - 1] << endl;
    return status;
  }

  int test_eval_poly()
  {
    Pruner<FT>::poly p;
    FT x, y;
    x      = .3;
    int ld = 3;
    p[0]   = 1;
    p[1]   = 1;
    p[2]   = 2;
    p[3]   = 3;

    y = pru.eval_poly(ld, p, x);
    cerr << "Testing eval_poly" << endl;
    cerr << y << endl;
    return abs(y - 1.561) > 1e-10;
  }
  int test_integrate_poly()
  {
    Pruner<FT>::poly p;
    int ld     = 3;
    int status = 0;
    p[0]       = 1.;
    p[1]       = 1.;
    p[2]       = 2.;
    p[3]       = 3.;
    pru.integrate_poly(ld, p);
    cerr << "Testing integrate_poly" << endl;
    for (int i = 0; i < 5; ++i)
    {
      cerr << p[i] << " ";
    }
    cerr << endl;
    status |= abs(p[0] - 0.0) > 1e-10;
    status |= abs(p[1] - 1.0) > 1e-10;
    status |= abs(p[2] - 1. / 2) > 1e-10;
    status |= abs(p[3] - 2. / 3) > 1e-10;
    status |= abs(p[4] - 3. / 4) > 1e-10;
    return status;
  }

  // To test relative volumes, we use the explicit formula
  // for step function given in [GNR], I_alpha(n/4,n/4)
  // computed using the following sage script
  //
  // alpha = XXX
  // n = 24
  // a = n/4
  // b = n/4
  // beta_dist = x**(a-1) * (1 - x)**(b-1)
  // c = integral(beta_dist, x, 0, alpha)
  // f = integral(beta_dist, x, 0, 1.)
  // print c/f

  // XXX = .3 ---> 0.07822479096
  // XXX = .5 ---> 0.5
  // XXX = .7 ---> 0.921775209040006

  int test_relative_volume()
  {
    vector<double> pr;
    pr.resize(Nbis);
    cerr << "Testing relative volume" << endl;
    int status = 0;

    double proba, error;

    for (int i = 0; i < Nbis / 2; ++i)
    {
      pr[i]              = 1;
      pr[i + (Nbis / 2)] = .3;
    }

    proba = fplll::svp_probability<FP_NR<double>>(pr).get_d();
    error = std::abs(1 - proba / 0.07822479096);
    cerr << proba << " relative error " << error << endl;
    status |= error > .05;

    for (int i = 0; i < Nbis / 2; ++i)
    {
      pr[i]              = 1;
      pr[i + (Nbis / 2)] = .5;
    }

    proba = fplll::svp_probability<FP_NR<double>>(pr).get_d();
    error = std::abs(1 - proba / 0.5);
    cerr << proba << " relative error " << error << endl;
    status |= error > .05;

    for (int i = 0; i < Nbis / 2; ++i)
    {
      pr[i]              = 1;
      pr[i + (Nbis / 2)] = .7;
    }

    proba = fplll::svp_probability<FP_NR<double>>(pr).get_d();
    error = std::abs(1 - proba / 0.92177520904);
    cerr << proba << " relative error " << error << endl;
    status |= error > .05;

    return status;
  }
};

void set_up_gso_norms(vector<double> &gso_sq_norms)
{
  for (int i = 0; i < N; ++i)
  {
    // A bit more than LLL reduced
    gso_sq_norms.emplace_back(pow(1.06, -i));
  }
}

template <class FT> int test_prepruned()
{
  int status = 0;
  cerr << "Checking Pre-pruned" << endl;
  Pruner<FT> pru;
  vector<double> gso_sq_norms;
  set_up_gso_norms(gso_sq_norms);
  pru.load_basis_shape(gso_sq_norms);

  vector<double> pr = {1,        1,        1,        1,        1,        1,        1,
                       1,        1,        1,        1,        1,        1,        1,
                       0.937197, 0.937197, 0.871731, 0.871731, 0.814304, 0.814304, 0.762232,
                       0.762232, 0.713898, 0.713898, 0.668279, 0.668279, 0.624701, 0.624701,
                       0.58271,  0.58271,  0.541994, 0.541994, 0.502342, 0.502342, 0.463617,
                       0.463617, 0.425747, 0.425747, 0.388723, 0.388723, 0.35262,  0.35262,
                       0.317642, 0.317642, 0.284261, 0.284261, 0.254584, 0.254584, 0.254584,
                       0.254584, 0.254584, 0.254584, 0.111895, 0.111895, 0.111895, 0.111895};

  pru.enumeration_radius = .85;
  double cost            = pru.single_enum_cost(pr);
  cerr << "Cost per enum " << cost << endl;
  status |= (abs(1 - cost / 1.7984e+07) > .01);
  double proba = pru.measure_metric(pr);
  cerr << "success proba " << proba << endl;
  status |= (abs(1 - proba / .506) > .01);
  return status;
}

template <class FT> int test_unpruned()
{
  int status = 0;
  cerr << "Checking Un-pruned" << endl;
  Pruner<FT> pru;
  vector<double> gso_sq_norms;
  set_up_gso_norms(gso_sq_norms);
  pru.load_basis_shape(gso_sq_norms);

  vector<double> pr;
  for (int i = 0; i < N; ++i)
  {
    pr.emplace_back(1.);
  }
  pru.enumeration_radius = .85;
  double cost            = pru.single_enum_cost(pr);
  cerr << "Cost per enum " << cost << endl;

  status       = (abs(1 - cost / 3.20e+10) > .02);
  double proba = pru.measure_metric(pr);
  cerr << "success proba " << proba << endl;
  status |= (abs(1 - proba) > .02);

  vector<vector<double>> gso_sq_norms_vec;
  vector<double> v1, v2, v3;
  v1 = gso_sq_norms;
  v2 = gso_sq_norms;
  v3 = gso_sq_norms;

  gso_sq_norms_vec.emplace_back(v1);
  gso_sq_norms_vec.emplace_back(v2);
  for (int i = 0; i < N; ++i)
  {
    v3[i] *= 20;
  }

  gso_sq_norms_vec.emplace_back(v3);
  pru.load_basis_shapes(gso_sq_norms_vec);

  cerr << "Repeating same checks with 3 bases" << endl;

  pru.enumeration_radius = .85;
  cost                   = pru.single_enum_cost(pr);
  cerr << "Cost per enum " << cost << endl;

  status = (abs(1 - 3. / 2. * cost / 3.20e+10) > .02);
  proba  = pru.measure_metric(pr);
  cerr << "success proba " << proba << endl;
  status |= (abs(1 - proba) > .02);
  return status;
}

template <class FT> int test_auto_prune(size_t n)
{
  int status = 0;
  IntMatrix A(2 * n, 2 * n);
  A.gen_ntrulike(30);
  IntMatrix U;
  MatGSO<Z_NR<mpz_t>, FP_NR<double>> M(A, U, U, GSO_DEFAULT);
  LLLReduction<Z_NR<mpz_t>, FP_NR<double>> lll_obj =
      LLLReduction<Z_NR<mpz_t>, FP_NR<double>>(M, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
  lll_obj.lll();
  FP_NR<double> radius;
  // NOTE: because NTRUlike lattice has a verri short vector 1111..
  // which is sometimes found by LLL, the pruner is only ran on dimension 1...2n-1.
  M.get_r(radius, 2, 2);
  vector<double> r;
  for (size_t i = 2; i < 2 * n; ++i)
  {
    FP_NR<double> x;
    M.get_r(x, i, i);
    r.push_back(x.get_d());
  }

  Pruning pruning;
  cerr << "Testing auto_prune " << endl;
  cerr << "RAD " << radius << endl;
  double radius_d = radius.get_d();

  double overhead = 1.0e6 * n * n;

  cerr << endl << "Gradient " << endl;
  prune<FT>(pruning, radius_d, overhead, 0.3, r, PRUNER_METHOD_GRADIENT,
            PRUNER_METRIC_PROBABILITY_OF_SHORTEST, true);
  status |= !(pruning.expectation <= 1.0);
  cerr << "Probability " << pruning.expectation << endl;
  status |= !(pruning.expectation > 0.0);
  status |= !(pruning.radius_factor >= 1.0);
  status |= !(pruning.coefficients[0] == 1.0);

  cerr << endl << "Reprune Gradient " << endl;
  prune<FT>(pruning, radius_d, overhead, 0.01, r, PRUNER_METHOD_GRADIENT,
            PRUNER_METRIC_PROBABILITY_OF_SHORTEST, false);
  status |= !(pruning.expectation <= 1.0);
  cerr << "Probability " << pruning.expectation << endl;
  status |= !(pruning.expectation > 0.0);
  status |= !(pruning.radius_factor >= 1.0);
  status |= !(pruning.coefficients[0] == 1.0);

  cerr << endl << "NelderMead " << endl;
  prune<FT>(pruning, radius_d, overhead, 0.3, r, PRUNER_METHOD_NM,
            PRUNER_METRIC_PROBABILITY_OF_SHORTEST, true);
  status |= !(pruning.expectation <= 1.0);
  cerr << "Probability " << pruning.expectation << endl;
  status |= !(pruning.expectation > 0.0);
  status |= !(pruning.radius_factor >= 1.0);
  status |= !(pruning.coefficients[0] == 1.0);

  cerr << endl << "Reprune NelderMead " << endl;
  prune<FT>(pruning, radius_d, overhead, 0.01, r, PRUNER_METHOD_GRADIENT,
            PRUNER_METRIC_PROBABILITY_OF_SHORTEST, false);
  status |= !(pruning.expectation <= 1.0);
  cerr << "Probability " << pruning.expectation << endl;
  status |= !(pruning.expectation > 0.0);
  status |= !(pruning.radius_factor >= 1.0);
  status |= !(pruning.coefficients[0] == 1.0);

  cerr << endl << "Hybrid " << endl;
  prune<FT>(pruning, radius_d, overhead, 0.3, r, PRUNER_METHOD_HYBRID,
            PRUNER_METRIC_PROBABILITY_OF_SHORTEST, true);
  status |= !(pruning.expectation <= 1.0);
  cerr << "Probability " << pruning.expectation << endl;
  status |= !(pruning.expectation > 0.0);
  status |= !(pruning.radius_factor >= 1.0);
  status |= !(pruning.coefficients[0] == 1.0);

  cerr << endl << "Reprune Hybrid " << endl;
  prune<FT>(pruning, radius_d, overhead, 0.01, r, PRUNER_METHOD_GRADIENT,
            PRUNER_METRIC_PROBABILITY_OF_SHORTEST, false);
  status |= !(pruning.expectation <= 1.0);
  cerr << "Probability " << pruning.expectation << endl;
  status |= !(pruning.expectation > 0.0);
  status |= !(pruning.radius_factor >= 1.0);
  status |= !(pruning.coefficients[0] == 1.0);

  cerr << endl << "Reprune Hybrid " << endl;
  prune<FT>(pruning, radius_d, overhead, 0.3, r, PRUNER_METHOD_GRADIENT,
            PRUNER_METRIC_PROBABILITY_OF_SHORTEST, false);
  status |= !(pruning.expectation <= 1.0);
  cerr << "Probability " << pruning.expectation << endl;
  status |= !(pruning.expectation > 0.0);
  status |= !(pruning.radius_factor >= 1.0);
  status |= !(pruning.coefficients[0] == 1.0);

  radius_d *= 2;
  cerr << endl << "Greedy " << endl;
  prune<FT>(pruning, radius_d, overhead, 20, r, PRUNER_METHOD_GREEDY,
            PRUNER_METRIC_EXPECTED_SOLUTIONS, true);
  status |= !(pruning.expectation > 1.0);
  cerr << "Probability " << pruning.expectation << endl;
  cerr << "Radius before/after " << 2 * radius.get_d() << "/" << radius_d << endl;
  status |= !(pruning.expectation > 0.0);
  status |= !(pruning.expectation < 100.0);
  status |= !(pruning.radius_factor >= 1.0);
  status |= !(pruning.coefficients[0] == 1.0);

  cerr << "NODES PER LEVEL" << endl;
  for (size_t i = 0; i < 2 * n - 1; ++i)
  {
    cerr << pruning.detailed_cost[i] << endl;
  }

  return status;
}

int main(int argc, char *argv[])
{

  int status = 0;
  status |= test_unpruned<FP_NR<double>>();
  status |= test_unpruned<FP_NR<long double>>();
#ifdef FPLLL_WITH_QD
  status |= test_unpruned<FP_NR<dd_real>>();
#endif
  status |= test_unpruned<FP_NR<mpfr_t>>();

  status |= test_prepruned<FP_NR<double>>();
  status |= test_prepruned<FP_NR<long double>>();
  status |= test_prepruned<FP_NR<mpfr_t>>();

  Pruner<FP_NR<long double>>::TestPruner tp;
  status |= tp.test_enforce();
  status |= tp.test_eval_poly();
  status |= tp.test_integrate_poly();
  status |= tp.test_relative_volume();

#ifdef FPLLL_WITH_QD
  Pruner<FP_NR<dd_real>>::TestPruner tp2;
  status |= tp2.test_enforce();
  status |= tp2.test_eval_poly();
  status |= tp2.test_integrate_poly();
  status |= tp2.test_relative_volume();
#endif

  status |= test_auto_prune<FP_NR<double>>(20);
  status |= test_auto_prune<FP_NR<double>>(30);

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
