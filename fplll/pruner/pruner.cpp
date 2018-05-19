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

template class Pruner<FP_NR<mpfr_t>>;

template void prune<FP_NR<mpfr_t>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template double prune_cost<FP_NR<mpfr_t>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<mpfr_t>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<mpfr_t> svp_probability<FP_NR<mpfr_t>>(const PruningParams &pruning);
template FP_NR<mpfr_t> svp_probability<FP_NR<mpfr_t>>(const vector<double> &pr);

// LD

#ifdef FPLLL_WITH_LONG_DOUBLE

template class Pruner<FP_NR<long double>>;
template void prune<FP_NR<long double>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template double prune_cost<FP_NR<long double>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<long double>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<long double> svp_probability<FP_NR<long double>>(const PruningParams &pruning);
template FP_NR<long double> svp_probability<FP_NR<long double>>(const vector<double> &pr);

#endif

#ifdef FPLLL_WITH_QD

template class Pruner<FP_NR<dd_real>>;
template void prune<FP_NR<dd_real>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template double prune_cost<FP_NR<dd_real>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<dd_real>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<dd_real> svp_probability<FP_NR<dd_real>>(const PruningParams &pruning);
template FP_NR<dd_real> svp_probability<FP_NR<dd_real>>(const vector<double> &pr);

template class Pruner<FP_NR<qd_real>>;
template void prune<FP_NR<qd_real>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<qd_real>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<qd_real> svp_probability<FP_NR<qd_real>>(const PruningParams &pruning);
template FP_NR<qd_real> svp_probability<FP_NR<qd_real>>(const vector<double> &pr);

#endif

#ifdef FPLLL_WITH_DPE

template class Pruner<FP_NR<dpe_t>>;
template void prune<FP_NR<dpe_t>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template double prune_cost<FP_NR<dpe_t>>(PruningParams &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<dpe_t>>(PruningParams &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<dpe_t> svp_probability<FP_NR<dpe_t>>(const PruningParams &pruning);
template FP_NR<dpe_t> svp_probability<FP_NR<dpe_t>>(const vector<double> &pr);

#endif
/* clang-format on */

FPLLL_END_NAMESPACE
