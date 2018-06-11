/* Copyright (C) 2015-2016 Martin Albrecht, Leo Ducas.

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

#include "pruner.h"
#include "ballvol.const"
#include "factorial.const"
#include "fplll.h"

// add components
#include "pruner_cost.cpp"
#include "pruner_optimize.cpp"
#include "pruner_optimize_tc.cpp"
#include "pruner_optimize_tp.cpp"
#include "pruner_prob.cpp"
#include "pruner_util.cpp"

FPLLL_BEGIN_NAMESPACE

// call pruner (note the float type is determined now)
template <class FT> int run_pruner_f(ZZ_mat<mpz_t> &b, const PruningParams &param, int sel_ft)
{
  int gso_flags = 0;
  if (b.get_rows() == 0 || b.get_cols() == 0)
    return RED_SUCCESS;
  if (sel_ft == FT_DOUBLE || sel_ft == FT_LONG_DOUBLE)
    gso_flags |= GSO_ROW_EXPO;

  // some checks
  int start = param.prune_start;
  int end   = param.prune_end;
  if (start < 0 || start >= b.get_rows() - 1)
    start = 0;
  if (end <= start || end >= b.get_rows())
    end                  = b.get_rows();
  double prune_pre_nodes = param.prune_pre_nodes;
  double prune_min_prob  = param.prune_min_prob;
  if (prune_pre_nodes <= 1)
    prune_pre_nodes = 1;
  int block_size    = end - start;

  PruningParams pruning;
  vector<double> r;
  FT root_det, max_dist;
  long max_dist_expo;

  // we check if we can convert the basis to long integers for
  // performance
  ZZ_mat<long> bl;
  if (convert<long, mpz_t>(bl, b, 10))
  {
    ZZ_mat<long> empty_mat;
    MatGSO<Z_NR<long>, FT> m_gso(bl, empty_mat, empty_mat, gso_flags);
    m_gso.update_gso();
    max_dist = m_gso.get_r_exp(start, start, max_dist_expo);
    root_det = m_gso.get_root_det(start, end);
    for (int i = start; i < end; ++i)
    {
      FT x;
      m_gso.get_r(x, i, i);
      r.push_back(x.get_d());
    }
  }
  else
  {
    ZZ_mat<mpz_t> empty_mat;
    MatGSO<Z_NR<mpz_t>, FT> m_gso(b, empty_mat, empty_mat, gso_flags);
    m_gso.update_gso();
    max_dist = m_gso.get_r_exp(start, start, max_dist_expo);
    root_det = m_gso.get_root_det(start, end);
    for (int i = start; i < end; ++i)
    {
      FT x;
      m_gso.get_r(x, i, i);
      r.push_back(x.get_d());
    }
  }

  adjust_radius_to_gh_bound(max_dist, max_dist_expo, block_size, root_det, param.gh_factor);
  double radius_d = max_dist.get_d() * pow(2, max_dist_expo);

  cerr << "# Start Pruning" << endl;
  cerr << "# enumeration Radius: " << radius_d << endl;
  cerr << "# preprocessing (num. nodes): " << prune_pre_nodes << endl;
  cerr << "# targeted min. prob: " << prune_min_prob << endl;
  cerr << "# input GSO: " << r << endl;
  prune<FT>(pruning, radius_d, prune_pre_nodes, r, prune_min_prob, PRUNER_METRIC_EXPECTED_SOLUTIONS,
            PRUNER_ZEALOUS | PRUNER_OPTIMIZE_FULL);
  cerr << "# optimized pruning coeff: " << endl << pruning.coefficients << endl;
  double cost = 0.;
  //  cerr << "# cost per level" << endl;
  for (int i = 0; i < block_size; ++i)
  {
    // cerr << pruning.detailed_cost[i] << " ";
    cost += pruning.detailed_cost[i];
  }
  cerr << "# single_enum_cost   = " << cost << endl;
  cerr << "#       succ. prob   = " << pruning.expectation << endl;
  cerr << "# repeated_enum_cost = " << cost / pruning.expectation << endl;
  return 0;
}

// interface function called from main.cpp
int run_pruner(ZZ_mat<mpz_t> &B, const PruningParams &param, FloatType float_type, int precision)
{
  // FPLLL_CHECK(B, "B == NULL in run_pruner()");
  FloatType sel_ft = (float_type != FT_DEFAULT) ? float_type : FT_DOUBLE;
  FPLLL_CHECK(!(sel_ft == FT_MPFR && precision == 0),
              "Missing precision for run_pruner() with floating point type mpfr");

  /* run pruner with float_type */
  int status;
  if (sel_ft == FT_DOUBLE)
  {
    status = run_pruner_f<FP_NR<double>>(B, param, sel_ft);
  }
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (sel_ft == FT_LONG_DOUBLE)
  {
    status = run_pruner_f<FP_NR<long double>>(B, param, sel_ft);
  }
#endif
#ifdef FPLLL_WITH_DPE
  else if (sel_ft == FT_DPE)
  {
    status = run_pruner_f<FP_NR<dpe_t>>(B, param, sel_ft);
  }
#endif
#ifdef FPLLL_WITH_QD
  else if (sel_ft == FT_DD)
  {
    status = run_pruner_f<FP_NR<dd_real>>(B, param, sel_ft);
  }
  else if (sel_ft == FT_QD)
  {
    status = run_pruner_f<FP_NR<qd_real>>(B, param, sel_ft);
  }
#endif
  else if (sel_ft == FT_MPFR)
  {
    int old_prec = FP_NR<mpfr_t>::set_prec(precision);
    status       = run_pruner_f<FP_NR<mpfr_t>>(B, param, sel_ft);
    FP_NR<mpfr_t>::set_prec(old_prec);
  }
  else
  {
    if (0 <= sel_ft && sel_ft <= FT_MPFR)
    {
      FPLLL_ABORT("Compiled without support for run_pruner() with " << FLOAT_TYPE_STR[sel_ft]);
    }
    else
    {
      FPLLL_ABORT("Floating point type " << sel_ft << "not supported in run_pruner()");
    }
  }
  return status;
}

template <class FT> FT svp_probability(const PruningParams &pruning)
{
  Pruner<FT> pru(pruning.coefficients.size());
  return pru.measure_metric(pruning.coefficients);
}

template <class FT> FT svp_probability(const vector<double> &pr)
{
  Pruner<FT> pru(pr.size());
  return pru.measure_metric(pr);
}

template <class FT> void Pruner<FT>::set_tabulated_consts()
{
  if (tabulated_values_imported)
    return;
  for (int i = 0; i < PRUNER_MAX_N; ++i)
  {
    tabulated_factorial[i] = pre_factorial[i];
    tabulated_ball_vol[i]  = pre_ball_vol[i];
  }
  tabulated_values_imported = 1;
  return;
}

template <class FT>
void prune(/*output*/ PruningParams &pruning,
           /*inputs*/ const double enumeration_radius, const double preproc_cost,
           const vector<double> &gso_r, const double target, const PrunerMetric metric,
           const int flags)
{
  Pruner<FT> pruner(enumeration_radius, preproc_cost, gso_r, target, metric, flags);
  pruner.optimize_coefficients(pruning.coefficients);
  pruner.single_enum_cost(pruning.coefficients, &(pruning.detailed_cost));
  pruning.gh_factor   = enumeration_radius / pruner.gaussian_heuristic().get_d();
  pruning.metric      = metric;
  pruning.expectation = pruner.measure_metric(pruning.coefficients);
}

template <class FT>
void prune(/*output*/ PruningParams &pruning,
           /*inputs*/ double enumeration_radius, const double preproc_cost,
           const vector<vector<double>> &gso_rs, const double target, const PrunerMetric metric,
           const int flags)
{
  Pruner<FT> pruner(enumeration_radius, preproc_cost, gso_rs, target, metric, flags);
  pruner.optimize_coefficients(pruning.coefficients);
  pruner.single_enum_cost(pruning.coefficients, &(pruning.detailed_cost));
  pruning.gh_factor   = enumeration_radius / pruner.gaussian_heuristic().get_d();
  pruning.metric      = metric;
  pruning.expectation = pruner.measure_metric(pruning.coefficients);
}

template <class FT>
double prune_cost(/*output*/ PruningParams &pruning,
                  /*inputs*/ const double enumeration_radius, const double preproc_cost,
                  const vector<double> &gso_r, const double target, const PrunerMetric metric,
                  const int flags)
{
  Pruner<FT> pruner(enumeration_radius, preproc_cost, gso_r, target, metric, flags);
  return pruner.repeated_enum_cost(pruning.coefficients);
}

/** instantiate functions **/
/* clang-format off */


// DOUBLE
template class Pruner<FP_NR<double>>;
template void prune<FP_NR<double>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template double prune_cost<FP_NR<double>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<double>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<double> svp_probability<FP_NR<double>>(const PruningParams &pruning);
template FP_NR<double> svp_probability<FP_NR<double>>(const vector<double> &pr);
template int run_pruner_f<FP_NR<double>> (ZZ_mat<mpz_t> &b, const PruningParams &param, int sel_ft);


// MPFR
template class Pruner<FP_NR<mpfr_t>>;
template void prune<FP_NR<mpfr_t>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template double prune_cost<FP_NR<mpfr_t>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<mpfr_t>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<mpfr_t> svp_probability<FP_NR<mpfr_t>>(const PruningParams &pruning);
template FP_NR<mpfr_t> svp_probability<FP_NR<mpfr_t>>(const vector<double> &pr);
template int run_pruner_f<FP_NR<mpfr_t>> (ZZ_mat<mpz_t> &b, const PruningParams &param, int sel_ft);


// LD
#ifdef FPLLL_WITH_LONG_DOUBLE

template class Pruner<FP_NR<long double>>;
template void prune<FP_NR<long double>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template double prune_cost<FP_NR<long double>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<long double>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<long double> svp_probability<FP_NR<long double>>(const PruningParams &pruning);
template FP_NR<long double> svp_probability<FP_NR<long double>>(const vector<double> &pr);
template int run_pruner_f<FP_NR<long double>> (ZZ_mat<mpz_t> &b, const PruningParams &param, int sel_ft);

#endif


#ifdef FPLLL_WITH_QD

// DD
template class Pruner<FP_NR<dd_real>>;
template void prune<FP_NR<dd_real>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template double prune_cost<FP_NR<dd_real>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<dd_real>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<dd_real> svp_probability<FP_NR<dd_real>>(const PruningParams &pruning);
template FP_NR<dd_real> svp_probability<FP_NR<dd_real>>(const vector<double> &pr);
template int run_pruner_f<FP_NR<dd_real>> (ZZ_mat<mpz_t> &b, const PruningParams &param, int sel_ft);


// QD
template class Pruner<FP_NR<qd_real>>;
template void prune<FP_NR<qd_real>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<qd_real>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<qd_real> svp_probability<FP_NR<qd_real>>(const PruningParams &pruning);
template FP_NR<qd_real> svp_probability<FP_NR<qd_real>>(const vector<double> &pr);
template int run_pruner_f<FP_NR<qd_real>> (ZZ_mat<mpz_t> &b, const PruningParams &param, int sel_ft);

#endif


#ifdef FPLLL_WITH_DPE
// DPE
template class Pruner<FP_NR<dpe_t>>;
template void prune<FP_NR<dpe_t>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template double prune_cost<FP_NR<dpe_t>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<dpe_t>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<dpe_t> svp_probability<FP_NR<dpe_t>>(const PruningParams &pruning);
template FP_NR<dpe_t> svp_probability<FP_NR<dpe_t>>(const vector<double> &pr);
template int run_pruner_f<FP_NR<dpe_t>> (ZZ_mat<mpz_t> &b, const PruningParams &param, int sel_ft);

#endif
/* clang-format on */

FPLLL_END_NAMESPACE
