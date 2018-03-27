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
#define D (N / 2)
#define Nbis 24
#define Dbis 12

int last_status = 0;
void print_status(int status)
{
  if (status > last_status)
  {
    last_status = status;
    cerr << endl << "FAILURES: " << status << endl;
  }
}

template <class FT> class Pruner<FT>::TestPruner
{
public:
  int n;
  int d;
  Pruner<FT> pru;
  TestPruner(int n) : n(n), pru(n) { d = n / 2; }

  int test_enforce()
  {
    vector<double> rs;
    int status = 0;
    for (int i = 0; i < n; ++i)
    {
      rs.emplace_back(1.);
    }
    pru.load_basis_shape(rs);
    Pruner<FT>::evec b(n);
    for (int i = 0; i < d; ++i)
    {
      b[i] = .3 + (1. * ((i * i * i) % d)) / (2 * d);
    }
    int j = d / 2;
    b[j]  = .5;

    FT old_bj = b[j];
    pru.min_pruning_coefficients.resize(d);
    fill(pru.min_pruning_coefficients.begin(), pru.min_pruning_coefficients.end(), .1);
    pru.enforce(b, j);

    status += !((b[j] >= old_bj) && (b[j] <= old_bj));
    print_status(status);
    status += !((b[d - 1] >= 1.) && (b[d - 1] <= 1.));
    print_status(status);
    cerr << "Checking enforced bounds: Increasing, 1-terminated, keeping b[13] = .5" << endl;
    for (int i = 0; i < d - 1; ++i)
    {
      status += !(b[i + 1] >= b[i]);
      print_status(status);
      cerr << b[i] << " , ";
    }
    cerr << b[d - 1] << endl;
    return status;
  }

  int test_eval_poly()
  {
    Pruner<FT>::poly p(4);
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
    Pruner<FT>::poly p(5);
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
    status += abs(p[0] - 0.0) > 1e-10;
    print_status(status);
    status += abs(p[1] - 1.0) > 1e-10;
    print_status(status);
    status += abs(p[2] - 1. / 2) > 1e-10;
    print_status(status);
    status += abs(p[3] - 2. / 3) > 1e-10;
    print_status(status);
    status += abs(p[4] - 3. / 4) > 1e-10;
    print_status(status);
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
    vector<double> pr(n);
    cerr << "Testing relative volume" << endl;
    int status = 0;

    double proba, error;

    for (int i = 0; i < n / 2; ++i)
    {
      pr[i]           = 1;
      pr[i + (n / 2)] = .3;
    }

    proba = fplll::svp_probability<FP_NR<double>>(pr).get_d();
    error = std::abs(1 - proba / 0.07822479096);
    cerr << proba << " relative error " << error << endl;
    status += error > .05;
    print_status(status);

    for (int i = 0; i < n / 2; ++i)
    {
      pr[i]           = 1;
      pr[i + (n / 2)] = .5;
    }

    proba = fplll::svp_probability<FP_NR<double>>(pr).get_d();
    error = std::abs(1 - proba / 0.5);
    cerr << proba << " relative error " << error << endl;
    status += error > .05;
    print_status(status);

    for (int i = 0; i < n / 2; ++i)
    {
      pr[i]           = 1;
      pr[i + (n / 2)] = .7;
    }

    proba = fplll::svp_probability<FP_NR<double>>(pr).get_d();
    error = std::abs(1 - proba / 0.92177520904);
    cerr << proba << " relative error " << error << endl;
    status += error > .05;
    print_status(status);

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
  cerr << endl << "Checking Pre-pruned" << endl;

  vector<double> gso_sq_norms;
  set_up_gso_norms(gso_sq_norms);

  PruningParams pruning;
  vector<double> pr = {1,        1,        1,        1,        1,        1,        1,
                       1,        1,        1,        1,        1,        1,        1,
                       0.937197, 0.937197, 0.871731, 0.871731, 0.814304, 0.814304, 0.762232,
                       0.762232, 0.713898, 0.713898, 0.668279, 0.668279, 0.624701, 0.624701,
                       0.58271,  0.58271,  0.541994, 0.541994, 0.502342, 0.502342, 0.463617,
                       0.463617, 0.425747, 0.425747, 0.388723, 0.388723, 0.35262,  0.35262,
                       0.317642, 0.317642, 0.284261, 0.284261, 0.254584, 0.254584, 0.254584,
                       0.254584, 0.254584, 0.254584, 0.2,      0.2,      0.2,      0.2};

  Pruner<FT> pru(.85, 20000., gso_sq_norms);
  double cost = pru.single_enum_cost(pr);
  cerr << "Cost per enum " << cost << endl;
  status += std::isnan(cost);
  print_status(status);
  status += (abs(1 - cost / 2.01206e+07) > .01);
  print_status(status);
  double proba = pru.measure_metric(pr);
  cerr << "success proba " << proba << endl;
  status += (abs(1 - proba / .556) > .01);
  print_status(status);
  status += std::isnan(proba);
  print_status(status);
  return status;
}

template <class FT> int test_unpruned()
{

  int status = 0;
  cerr << "Checking Un-pruned" << endl;
  vector<double> gso_sq_norms;
  set_up_gso_norms(gso_sq_norms);

  vector<double> pr;
  for (int i = 0; i < N; ++i)
  {
    pr.emplace_back(1.);
  }
  Pruner<FT> pru(.85, 20000., gso_sq_norms);

  double cost = pru.single_enum_cost(pr);
  cerr << "Cost per enum " << cost << endl;
  status += (abs(1 - cost / 3.20e+10) > .02);
  print_status(status);
  status += std::isnan(cost);
  print_status(status);
  double proba = pru.measure_metric(pr);
  cerr << "success proba " << proba << endl;
  status += (abs(1 - proba) > .02);
  print_status(status);
  status += std::isnan(proba);
  print_status(status);

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
  // 1 basis out of 3 is so large that it essentially induces no cost: new cost should be 2/3 of
  // the previous

  gso_sq_norms_vec.emplace_back(v3);

  cerr << "Repeating same checks with 3 bases" << endl;
  Pruner<FT> pru2(.85, 0., gso_sq_norms_vec);

  cost = pru2.single_enum_cost(pr);
  cerr << "Cost per enum " << cost << endl;

  status += (abs(1 - 3. / 2. * cost / 3.20e+10) >
             .02);  // check that the new cost is indeed 2/3 of the original cost
  print_status(status);
  proba = pru2.measure_metric(pr);
  cerr << "success proba " << proba << endl;
  status += (abs(1 - proba) > .02);
  print_status(status);
  return status;
  return 0;
}

template <class FT> int test_auto_prune(size_t n)
{
  int status = 0;
  double cost;
  ZZ_mat<mpz_t> A(2 * n, 2 * n);
  A.gen_qary(n, 30);
  ZZ_mat<mpz_t> U;
  MatGSO<Z_NR<mpz_t>, FP_NR<double>> M(A, U, U, GSO_DEFAULT);
  LLLReduction<Z_NR<mpz_t>, FP_NR<double>> lll_obj =
      LLLReduction<Z_NR<mpz_t>, FP_NR<double>>(M, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
  lll_obj.lll();
  FP_NR<double> radius;
  M.get_r(radius, 0, 0);
  vector<double> r;
  for (size_t i = 0; i < 2 * n; ++i)
  {
    FP_NR<double> x;
    M.get_r(x, i, i);
    r.push_back(x.get_d());
  }

  PruningParams pruning;
  cerr << "Testing auto_prune " << endl;
  double overhead = 1.0e6 * n * n;
  cerr << "Overhead " << overhead << endl;

  double radius_d = r[0] * .3;

  cerr << endl << "Greedy " << endl;
  prune<FT>(pruning, radius_d, overhead, r, 20, PRUNER_METRIC_EXPECTED_SOLUTIONS, 0);
  cerr << "Expected Solutions " << pruning.expectation << endl;
  cerr << "Radius " << radius_d << endl;
  cerr << "gh_factor " << pruning.gh_factor << endl;

  status += !(pruning.expectation > 0.0);
  print_status(status);
  status += !(pruning.gh_factor >= .05);
  print_status(status);
  status += !(pruning.gh_factor < 20.);
  print_status(status);

  status += !(pruning.coefficients[0] == 1.0);
  print_status(status);
  cost = 0.;
  cerr << "Predicted cost per Level" << endl;
  for (size_t i = 0; i < 2 * n; ++i)
  {
    cerr << pruning.detailed_cost[i] << "\t";
    cost += pruning.detailed_cost[i];
  }
  cerr << endl << "Predicted Total Cost " << cost << endl;

  cerr << endl << "Gradient " << endl;
  cerr << "radius " << radius_d << endl;
  prune<FT>(pruning, radius_d, overhead, r, 0.3, PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
            PRUNER_GRADIENT);
  status += !(pruning.expectation <= 1.001);
  print_status(status);
  cerr << "Probability " << pruning.expectation << endl;
  cost = 0.;
  cerr << "Predicted cost per Level" << endl;
  for (size_t i = 0; i < 2 * n; ++i)
  {
    cerr << pruning.detailed_cost[i] << "\t";
    cost += pruning.detailed_cost[i];
  }
  cerr << endl << "Predicted Total Cost " << cost << endl;

  status += !(pruning.expectation > 0.0);
  print_status(status);
  status += !(pruning.gh_factor >= .05);
  print_status(status);
  status += !(pruning.gh_factor < 20.);
  print_status(status);
  status += !(pruning.coefficients[0] == 1.0);
  print_status(status);

  cerr << endl << "Reprune Gradient " << endl;
  cerr << "radius " << radius_d << endl;
  prune<FT>(pruning, radius_d, overhead, r, 0.01, PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
            PRUNER_GRADIENT | PRUNER_START_FROM_INPUT);
  status += !(pruning.expectation <= 1.001);
  print_status(status);
  cerr << "Probability " << pruning.expectation << endl;
  cost = 0.;
  cerr << "Predicted cost per Level" << endl;
  for (size_t i = 0; i < 2 * n; ++i)
  {
    cerr << pruning.detailed_cost[i] << "\t";
    cost += pruning.detailed_cost[i];
  }
  cerr << endl << "Predicted Total Cost " << cost << endl;
  status += !(pruning.expectation > 0.0);
  print_status(status);
  status += !(pruning.gh_factor >= .05);
  print_status(status);
  status += !(pruning.gh_factor < 20.);
  print_status(status);
  status += !(pruning.coefficients[0] == 1.0);
  print_status(status);

  cerr << endl << "NelderMead " << endl;
  cerr << "radius " << radius_d << endl;
  prune<FT>(pruning, radius_d, overhead, r, 0.3, PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
            PRUNER_NELDER_MEAD);
  status += !(pruning.expectation <= 1.001);
  print_status(status);
  cerr << "Probability " << pruning.expectation << endl;
  cost = 0.;
  cerr << "Predicted cost per Level" << endl;
  for (size_t i = 0; i < 2 * n; ++i)
  {
    cerr << pruning.detailed_cost[i] << "\t";
    cost += pruning.detailed_cost[i];
  }
  cerr << endl << "Predicted Total Cost " << cost << endl;
  status += !(pruning.expectation > 0.0);
  print_status(status);
  status += !(pruning.gh_factor >= .05);
  print_status(status);
  status += !(pruning.gh_factor < 20.);
  print_status(status);
  status += !(pruning.coefficients[0] == 1.0);
  print_status(status);

  cerr << endl << "Reprune NelderMead " << endl;
  cerr << "radius " << radius_d << endl;
  prune<FT>(pruning, radius_d, overhead, r, 0.3, PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
            PRUNER_NELDER_MEAD | PRUNER_START_FROM_INPUT);
  status += !(pruning.expectation <= 1.001);
  print_status(status);
  cerr << "Probability " << pruning.expectation << endl;
  cost = 0.;
  cerr << "Predicted cost per Level" << endl;
  for (size_t i = 0; i < 2 * n; ++i)
  {
    cerr << pruning.detailed_cost[i] << "\t";
    cost += pruning.detailed_cost[i];
  }
  cerr << endl << "Predicted Total Cost " << cost << endl;
  status += !(pruning.expectation > 0.0);
  print_status(status);
  status += !(pruning.gh_factor >= .05);
  print_status(status);
  status += !(pruning.gh_factor < 20.);
  print_status(status);

  status += !(pruning.coefficients[0] == 1.0);
  print_status(status);

  cerr << endl << "Zealous (Gradient then NelderMead) " << endl;
  cerr << "radius " << radius_d << endl;
  prune<FT>(pruning, radius_d, overhead, r, 0.3, PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
            PRUNER_ZEALOUS);
  status += !(pruning.expectation <= 1.001);
  print_status(status);
  cerr << "Probability " << pruning.expectation << endl;
  cost = 0.;
  cerr << "Predicted cost per Level" << endl;
  for (size_t i = 0; i < 2 * n; ++i)
  {
    cerr << pruning.detailed_cost[i] << "\t";
    cost += pruning.detailed_cost[i];
  }
  cerr << endl << "Predicted Total Cost " << cost << endl;
  status += !(pruning.expectation > 0.0);
  print_status(status);
  status += !(pruning.gh_factor >= .05);
  print_status(status);
  status += !(pruning.gh_factor < 20.);
  print_status(status);

  status += !(pruning.coefficients[0] == 1.0);
  print_status(status);

  cerr << endl << "Reprune Zealous " << endl;
  cerr << "radius " << radius_d << endl;
  prune<FT>(pruning, radius_d, overhead, r, 0.3, PRUNER_METRIC_PROBABILITY_OF_SHORTEST,
            PRUNER_ZEALOUS | PRUNER_START_FROM_INPUT);
  status += !(pruning.expectation <= 1.001);
  print_status(status);
  cerr << "Probability " << pruning.expectation << endl;
  cost = 0.;
  cerr << "Predicted cost per Level" << endl;
  for (size_t i = 0; i < 2 * n; ++i)
  {
    cerr << pruning.detailed_cost[i] << "\t";
    cost += pruning.detailed_cost[i];
  }
  cerr << endl << "Predicted Total Cost " << cost << endl;
  status += !(pruning.expectation > 0.0);
  print_status(status);
  status += !(pruning.gh_factor >= .05);
  print_status(status);
  status += !(pruning.gh_factor < 20.);
  print_status(status);

  status += !(pruning.coefficients[0] == 1.0);
  print_status(status);

  return status;
}

int main()
{
  int status = 0;
#ifdef FPLLL_WITH_QD
  cerr << endl << "DD" << endl;
  status += test_unpruned<FP_NR<dd_real>>();
  print_status(status);

#endif

  cerr << endl << "d" << endl;
  status += test_unpruned<FP_NR<double>>();
  print_status(status);
#ifdef FPLLL_WITH_LONG_DOUBLE
  cerr << endl << "ld" << endl;
  status += test_unpruned<FP_NR<long double>>();
  print_status(status);
#endif
  cerr << endl << "MPRF" << endl;
  status += test_unpruned<FP_NR<mpfr_t>>();
  print_status(status);

  status += test_prepruned<FP_NR<double>>();
  print_status(status);
#ifdef FPLLL_WITH_LONG_DOUBLE
  status += test_prepruned<FP_NR<long double>>();
  print_status(status);
#endif
  status += test_prepruned<FP_NR<mpfr_t>>();
  print_status(status);

#ifdef FPLLL_WITH_LONG_DOUBLE
  Pruner<FP_NR<long double>>::TestPruner tp(Nbis);
  status += tp.test_enforce();
  print_status(status);
  status += tp.test_eval_poly();
  print_status(status);
  status += tp.test_integrate_poly();
  print_status(status);
  status += tp.test_relative_volume();
  print_status(status);
#endif

#ifdef FPLLL_WITH_QD
  Pruner<FP_NR<dd_real>>::TestPruner tp2(Nbis);
  status += tp2.test_enforce();
  print_status(status);
  status += tp2.test_eval_poly();
  print_status(status);
  status += tp2.test_integrate_poly();
  print_status(status);
  status += tp2.test_relative_volume();
  print_status(status);
#endif

#ifdef FPLLL_WITH_QD
  status += test_auto_prune<FP_NR<dd_real>>(20);
  print_status(status);
#endif

  status += test_auto_prune<FP_NR<double>>(20);
  print_status(status);
  status += test_auto_prune<FP_NR<double>>(30);
  print_status(status);

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
