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
#include "gso.h"

FPLLL_BEGIN_NAMESPACE

// PRIVATE METHODS

template <class FT> FT svp_probability(const Pruning &pruning)
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

template <class FT> void Pruner<FT>::optimize_coefficients(/*io*/ vector<double> &pr)
{
  evec b(d);
  if (flags & PRUNER_START_FROM_INPUT)
  {
    load_coefficients(b, pr);
  }
  if (!(flags & PRUNER_START_FROM_INPUT))
  {
    greedy(b);
  }

  if (flags & (PRUNER_GRADIENT | PRUNER_NELDER_MEAD))
  {
    preproc_cost *= .1;
    greedy(min_pruning_coefficients);
    preproc_cost *= 10;
  }

  if (flags & PRUNER_GRADIENT)
  {
    while (gradient_descent_step(b))
    {
    };
  };
  if (flags & PRUNER_NELDER_MEAD)
  {
    while (nelder_mead_step(b))
    {
    };
  };
  save_coefficients(pr, b);
}

template <class FT>
void Pruner<FT>::load_basis_shape(const vector<double> &gso_r, bool reset_normalization)
{
  shape_loaded = true;
  FT logvol, tmp;
  logvol = 0.0;
  r.resize(n);
  ipv.resize(n);
  for (int i = 0; i < n; ++i)
  {
    r[i] = gso_r[n - 1 - i];
    logvol += log(r[i]);
  }

  if (reset_normalization)
  {
    normalization_factor = exp(logvol / ((float)(-n)));
    normalized_radius    = sqrt(enumeration_radius * normalization_factor);
  }

  for (int i = 0; i < n; ++i)
  {
    r[i] *= normalization_factor;
  }
  tmp = 1.;
  for (int i = 0; i < 2 * d; ++i)
  {
    tmp *= sqrt(r[i]);
    ipv[i] = 1.0 / tmp;
  }
}

template <class FT> void Pruner<FT>::load_basis_shapes(const vector<vector<double>> &gso_rs)
{
  n = gso_rs[0].size();
  vec sum_ipv(n);

  for (int i = 0; i < n; ++i)
  {
    sum_ipv[i] = 0.;
  }
  int count = gso_rs.size();
  for (int k = 0; k < count; ++k)
  {
    if (gso_rs[k].size() != static_cast<size_t>(n))
    {
      throw std::runtime_error("loading several bases with different dimensions");
    }
    load_basis_shape(gso_rs[k], (k == 0));
    for (int i = 0; i < n; ++i)
    {
      sum_ipv[i] += ipv[i];
    }
  }
  for (int i = 0; i < n; ++i)
  {
    ipv[i] = sum_ipv[i] / (1.0 * count);
  }
}

template <class FT> FT Pruner<FT>::gaussian_heuristic()
{
  return exp(2. * log(tabulated_ball_vol[n]) / ((float)-n)) / normalization_factor;
}

template <class FT>
void Pruner<FT>::load_coefficients(/*o*/ evec &b, /*i*/ const vector<double> &pr)
{
  for (int i = 0; i < d; ++i)
  {
    b[i] = pr[n - 1 - 2 * i];
  }
  if (enforce(b))
  {
    throw std::runtime_error(
        "Ill formed pruning coefficients (must be decreasing, starting with two 1.0)");
  }
}

template <class FT>
void Pruner<FT>::save_coefficients(/*o*/ vector<double> &pr, /*i*/ const evec &b)
{
  pr.resize(n);
  for (int i = 0; i < d; ++i)
  {
    pr[n - 1 - 2 * i] = b[i].get_d();
    pr[n - 2 - 2 * i] = b[i].get_d();
  }
  pr[0] = 1.;
}

template <class FT> inline bool Pruner<FT>::enforce(/*io*/ evec &b, /*opt i*/ const int j)
{
  bool status = false;
  if ((b[d - 1] < .999) & (j != d - 1))
  {
    status   = 1;
    b[d - 1] = 1.;
  }
  for (int i = 0; i < d; ++i)
  {
    status |= (b[i] > 1.0001);
    b[i] = b[i] > 1 ? 1. : b[i];

    if (b[i] <= min_pruning_coefficients[i])
      b[i] = min_pruning_coefficients[i];
  }
  for (int i = j; i < d - 1; ++i)
  {
    if (b[i + 1] < b[i])
    {
      status |= (b[i + 1] + .001 < b[i]);
      b[i + 1] = b[i];
    }
  }
  for (int i = j - 1; i >= 0; --i)
  {
    if (b[i + 1] < b[i])
    {
      status |= (b[i + 1] + .001 < b[i]);
      b[i] = b[i + 1];
    }
  }
  return status;
}

template <class FT> inline FT Pruner<FT>::eval_poly(const int ld, /*i*/ const poly &p, const FT x)
{
  FT acc;
  acc = 0.0;
  for (int i = ld; i >= 0; --i)
  {
    acc = acc * x;
    acc = acc + p[i];
  }
  return acc;
}

template <class FT> inline void Pruner<FT>::integrate_poly(const int ld, /*io*/ poly &p)
{
  for (int i = ld; i >= 0; --i)
  {
    FT tmp;
    tmp      = i + 1.;
    p[i + 1] = p[i] / tmp;
  }
  p[0] = 0.0;
}

template <class FT> inline FT Pruner<FT>::relative_volume(const int rd, /*i*/ const evec &b)
{
  poly P(rd + 1);
  P[0]   = 1;
  int ld = 0;
  for (int i = rd - 1; i >= 0; --i)
  {
    integrate_poly(ld, P);
    ld++;
    P[0] = -1.0 * eval_poly(ld, P, b[i] / b[rd - 1]);
  }
  FT res = P[0] * tabulated_factorial[rd];
  return (rd % 2) ? -res : res;
}

template <class FT>
inline FT Pruner<FT>::single_enum_cost(/*i*/ const evec &b, vector<double> *detailed_cost)
{
  if (!shape_loaded)
  {
    throw std::invalid_argument("No basis shape was loaded");
  }

  if (detailed_cost)
  {
    detailed_cost->resize(n);
  }
  vec rv(n);  // Relative volumes at each level

  for (int i = 0; i < d; ++i)
  {

    rv[2 * i + 1] = relative_volume(i + 1, b);
  }

  rv[0] = 1;
  for (int i = 1; i < d; ++i)
  {
    rv[2 * i] = sqrt(rv[2 * i - 1] * rv[2 * i + 1]);  // Interpolate even values
  }

  FT total;
  total                    = 0.0;
  FT normalized_radius_pow = normalized_radius;
  for (int i = 0; i < 2 * d; ++i)
  {
    FT tmp;

    tmp = normalized_radius_pow * rv[i] * tabulated_ball_vol[i + 1] *
          sqrt(pow_si(b[i / 2], 1 + i)) * ipv[i];
    tmp *= symmetry_factor;
    if (detailed_cost)
    {
      (*detailed_cost)[2 * d - (i + 1)] = tmp.get_d();
    }

    total += tmp;
    normalized_radius_pow *= normalized_radius;
  }
  return total;
}

template <class FT> inline FT Pruner<FT>::svp_probability(/*i*/ const evec &b)
{
  evec b_minus_db(d);
  FT dx = shell_ratio;
  for (int i = 0; i < d; ++i)
  {
    b_minus_db[i] = b[i] / (dx * dx);
    if (b_minus_db[i] > 1)
      b_minus_db[i] = 1;
  }

  FT vol  = relative_volume(d, b);
  FT dxn  = pow_si(dx, 2 * d);
  FT dvol = dxn * relative_volume(d, b_minus_db) - vol;
  return dvol / (dxn - 1.);
}

template <class FT> inline FT Pruner<FT>::expected_solutions(/*i*/ const evec &b)
{
  if (!shape_loaded)
  {
    throw std::invalid_argument("No basis shape was loaded");
  }

  int j  = d * 2 - 1;
  FT tmp = relative_volume((j + 1) / 2, b);
  tmp *= tabulated_ball_vol[j + 1];
  tmp *= pow_si(normalized_radius * sqrt(b[j / 2]), j + 1);
  tmp *= ipv[j];
  tmp *= symmetry_factor;

  return tmp;
}

template <class FT> inline FT Pruner<FT>::measure_metric(/*i*/ const evec &b)
{
  if (metric == PRUNER_METRIC_PROBABILITY_OF_SHORTEST)
  {
    return svp_probability(b);
  }
  else if (metric == PRUNER_METRIC_EXPECTED_SOLUTIONS)
  {
    return expected_solutions(b);
  }
  else
  {
    throw std::invalid_argument("Pruner was set to an unknown metric");
  }
}

template <class FT> inline FT Pruner<FT>::repeated_enum_cost(/*i*/ const evec &b)
{

  if (metric == PRUNER_METRIC_PROBABILITY_OF_SHORTEST)
  {
    FT probability = svp_probability(b);
    if (probability >= target)
      return single_enum_cost(b);

    FT trials = log(1.0 - target) / log(1.0 - probability);
    return single_enum_cost(b) * trials + preproc_cost * (trials - 1.0);
  }

  else if (metric == PRUNER_METRIC_EXPECTED_SOLUTIONS)
  {
    FT expected = expected_solutions(b);
    if (expected >= target)
      return single_enum_cost(b);

    FT trials = target / expected;
    if (trials < 1.)
      trials = 1;
    return single_enum_cost(b) * trials + preproc_cost * (trials - 1.0);
  }

  else
  {
    throw std::invalid_argument("Pruner was set to an unknown metric");
  }
}

template <class FT>
void Pruner<FT>::repeated_enum_cost_gradient(/*i*/ const evec &b, /*o*/ evec &res)
{
  evec b_plus_db(d);
  res[d - 1] = 0.0;  // Force null gradient on the last coordinate : don't touch this coeff
  for (int i = 0; i < d - 1; ++i)
  {
    b_plus_db = b;
    b_plus_db[i] *= (1.0 - epsilon);
    enforce(b_plus_db, i);
    FT X = repeated_enum_cost(b_plus_db);

    b_plus_db = b;
    b_plus_db[i] *= (1.0 + epsilon);
    enforce(b_plus_db, i);
    FT Y   = repeated_enum_cost(b_plus_db);
    res[i] = (log(X) - log(Y)) / epsilon;
  }
}

template <class FT> int Pruner<FT>::gradient_descent_step(/*io*/ evec &b)
{

  FT cf     = repeated_enum_cost(b);
  FT old_cf = cf;
  evec new_b(d);
  evec gradient(d);
  repeated_enum_cost_gradient(b, gradient);
  FT norm = 0.0;

  // normalize the gradient
  for (int i = 0; i < d; ++i)
  {
    norm += gradient[i] * gradient[i];
    new_b[i] = b[i];
  }

  norm /= (double)d;
  norm = sqrt(norm);
  if (norm <= 0.)
    return 0;

  for (int i = 0; i < d; ++i)
  {
    gradient[i] /= norm;
  }
  FT new_cf;

  FT step = min_step;
  int i;

  for (i = 0;; ++i)
  {
    for (int i = 0; i < d; ++i)
    {
      new_b[i] = new_b[i] + step * gradient[i];
    }

    enforce(new_b);
    new_cf = repeated_enum_cost(new_b);

    if (new_cf >= cf)
    {
      break;
    }
    b  = new_b;
    cf = new_cf;
    step *= step_factor;
  }
  if (cf > old_cf * min_cf_decrease)
  {
    return 0;
  }
  return i;
}

template <class FT> void Pruner<FT>::greedy(evec &b)
{
  // Do not call enforce in this function, as min_pruning_bounds may not have been set
  // Indeed, the min_pruning_bound should now based on greedy.
  if (!shape_loaded)
  {
    throw std::invalid_argument("No basis shape was loaded");
  }

  fill(min_pruning_coefficients.begin(), min_pruning_coefficients.end(), 0.);

  b.resize(d);
  fill(b.begin(), b.end(), 1.);
  evec new_b(d);
  FT nodes;

  for (int j = 1; j < 2 * d - 1; j += 2)
  {
    int i = j / 2;
    if (i > 1)
    {
      b[i] = b[i - 1] > .9 ? 1 : 1.1 * b[i - 1];
    }
    double goal_factor =
        1. / (3. * n) +
        4 * j * (n - j) / (n * n * n);  // Make the tree width as a parabola, with maximum at n/2
    nodes = 1. + 1e10 * preproc_cost;
    while ((nodes > goal_factor * preproc_cost) & (b[i] > .001))
    {
      b[i] *= .98;
      for (int k = 0; k < i; ++k)
      {
        b[k] = b[k] < b[i] ? b[k] : b[i];  // Enforcing decreasing by hand
      }
      nodes = relative_volume((j + 1) / 2, b);
      nodes *= tabulated_ball_vol[j + 1];
      nodes *= pow_si(normalized_radius * sqrt(b[i]), j + 1);
      nodes *= ipv[j];
      nodes *= symmetry_factor;
    }
  }
}

// Nelder-Mead method. Following the notation of
// https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method

#define ND_ALPHA 1
#define ND_GAMMA 2
#define ND_RHO 0.5
#define ND_SIGMA 0.5
#define ND_INIT_WIDTH 0.05

template <class FT> int Pruner<FT>::nelder_mead_step(/*io*/ evec &b)
{
  int l = d + 1;
  evec tmp_constructor(d);
  vector<evec> bs(l);  // The simplexe (d+1) vector of dim d
  FT *fs = new FT[l];  // Values of f at the simplex vertices

  for (int i = 0; i < l; ++i)  // Intialize the simplex
  {
    bs[i] = b;  // Start from b
    if (i < d)
    {
      bs[i][i] += (bs[i][i] < .5) ? ND_INIT_WIDTH : -ND_INIT_WIDTH;
    }
    enforce(bs[i]);
    fs[i] = repeated_enum_cost(bs[i]);  // initialize the value
  }

  FT init_cf = fs[l - 1];

  evec bo(d);  // centeroid
  FT fo;       // value at the centroid

  FT fs_maxi_last;  // value of the last centroid

  if (verbosity)
  {
    cerr << "  Starting nelder_mead cf = " << init_cf << " proba = " << svp_probability(b) << endl;
  }
  unsigned int counter = 0;
  int mini = 0, maxi = 0, maxi2 = 0;
  while (1)  // Main loop
  {
    mini = maxi = maxi2 = 0;
    for (int i = 0; i < d; ++i)
      bo[i]    = bs[0][i];
    ////////////////
    // step 1. and 2. : Order and centroid
    ////////////////
    for (int i = 1; i < l; ++i)  // determine min and max, and centroid
    {
      if (fs[i] < fs[mini])
        mini = i;
      if (fs[i] > fs[maxi])
        maxi = i;
      for (int j = 0; j < d; ++j)
        bo[j] += bs[i][j];
    }
    FT tmp;
    tmp = l;
    for (int i = 0; i < d; ++i)
      bo[i] /= tmp;  // Centroid calculated

    if (!counter)
      fs_maxi_last = fs[maxi];

    if (!maxi)
      maxi2++;
    for (int i = 1; i < l; ++i)  // determine min and max, and centroid
    {
      if ((fs[i] > fs[maxi2]) && (i != maxi))
        maxi2 = i;
    }

    if (enforce(bo))
    {
      throw std::runtime_error("Concavity says that should not happen.");
    }

    if (verbosity)
    {
      cerr << "  melder_mead step " << counter << "cf = " << fs[mini]
           << " proba = " << measure_metric(bs[mini]) << " cost = " << single_enum_cost(bs[mini])
           << endl;
      for (int i = 0; i < d; ++i)
      {
        cerr << ceil(bs[mini][i].get_d() * 1000) << " ";
      }
      cerr << endl;
    }

    ////////////////
    // Stopping condition (Not documented on wikipedia, improvising)
    ////////////////

    counter++;
    if (!(counter % l))  // Every l steps, we check progress and stop if none is done
    {
      if (fs[maxi] > fs_maxi_last * min_cf_decrease)
      {
        break;
      }
      fs_maxi_last = fs[maxi];
    }

    for (int i = 0; i < l; ++i)  // determine second best
    {
      if ((fs[i] > fs[maxi2]) && (i != maxi))
        maxi2 = i;
    }

    if (verbosity)
    {
      cerr << mini << " " << maxi2 << " " << maxi << endl;
      cerr << fs[mini] << " < " << fs[maxi2] << " < " << fs[maxi] << " | " << endl;
    }

    ////////////////
    // step 3. Reflection
    ////////////////

    evec br(d);  // reflected point
    FT fr;       // Value at the reflexion point
    for (int i = 0; i < d; ++i)
      br[i]    = bo[i] + ND_ALPHA * (bo[i] - bs[maxi][i]);
    enforce(br);
    fr = repeated_enum_cost(br);
    if (verbosity)
    {
      cerr << "fr " << fr << endl;
    }

    if ((fs[mini] <= fr) && (fr < fs[maxi2]))
    {
      bs[maxi] = br;
      fs[maxi] = fr;
      if (verbosity)
      {
        cerr << "    Reflection " << endl;
      }
      continue;  // Go to step 1.
    }

    ////////////////
    // step 4. Expansion
    ////////////////

    if (fr < fs[mini])
    {
      evec be(d);
      FT fe;
      for (int i = 0; i < d; ++i)
        be[i]    = bo[i] + ND_GAMMA * (br[i] - bo[i]);
      enforce(be);
      fe = repeated_enum_cost(be);
      if (verbosity)
      {
        cerr << "fe " << fe << endl;
      }
      if (fe < fr)
      {
        bs[maxi] = be;
        fs[maxi] = fe;
        if (verbosity)
        {
          cerr << "    Expansion A " << endl;
        }
        continue;  // Go to step 1.
      }
      else
      {
        bs[maxi] = br;
        fs[maxi] = fr;
        if (verbosity)
        {
          cerr << "    Expansion B " << endl;
        }
        continue;  // Go to step 1.
      }
    }

    ////////////////
    // step 5. Contraction
    ////////////////

    if (!(fr >= fs[maxi2]))  // Here, it is certain that fr >= fs[maxi2]
    {
      throw std::runtime_error("Something certain is false in Nelder-Mead.");
    }

    evec bc(d);
    FT fc;
    for (int i = 0; i < d; ++i)
      bc[i]    = bo[i] + ND_RHO * (bs[maxi][i] - bo[i]);
    enforce(bc);
    fc = repeated_enum_cost(bc);
    if (verbosity)
    {
      cerr << "fc " << fc << endl;
    }
    if (fc < fs[maxi])
    {
      bs[maxi] = bc;
      fs[maxi] = fc;
      if (verbosity)
      {
        cerr << "    Contraction " << endl;
      }
      continue;  // Go to step 1.
    }

    ////////////////
    // step 6. Shrink
    ////////////////
    if (verbosity)
    {
      cerr << "    Shrink " << endl;
    }
    for (int j = 0; j < l; ++j)
    {
      for (int i = 0; i < d; ++i)
      {
        bs[j][i] = bs[mini][i] + ND_SIGMA * (bs[j][i] - bs[mini][i]);
      }
      enforce(bs[j]);
      fs[j] = repeated_enum_cost(bs[j]);  // initialize the value
    }
  }

  b            = bs[mini];
  int improved = (init_cf * min_cf_decrease) > fs[mini];

  if (verbosity)
  {
    cerr << "Done nelder_mead, after " << counter << " steps" << endl;
    cerr << "Final cf = " << fs[mini] << " proba = " << measure_metric(b) << endl;
    if (improved)
    {
      cerr << "Progress has been made: init cf = " << init_cf << endl;
    }
    cerr << endl;
  }

  return improved;  // Has MN made any progress
}

template <class FT>
void prune(/*output*/ Pruning &pruning,
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
void prune(/*output*/ Pruning &pruning,
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

/** instantiate functions **/
/* clang-format off */

// DOUBLE

template class Pruner<FP_NR<double>>;

template void prune<FP_NR<double>>(Pruning &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<double>>(Pruning &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<double> svp_probability<FP_NR<double>>(const Pruning &pruning);
template FP_NR<double> svp_probability<FP_NR<double>>(const vector<double> &pr);

template class Pruner<FP_NR<mpfr_t>>;

template void prune<FP_NR<mpfr_t>>(Pruning &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<mpfr_t>>(Pruning &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<mpfr_t> svp_probability<FP_NR<mpfr_t>>(const Pruning &pruning);
template FP_NR<mpfr_t> svp_probability<FP_NR<mpfr_t>>(const vector<double> &pr);

// LD

#ifdef FPLLL_WITH_LONG_DOUBLE

template class Pruner<FP_NR<long double>>;
template void prune<FP_NR<long double>>(Pruning &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<long double>>(Pruning &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<long double> svp_probability<FP_NR<long double>>(const Pruning &pruning);
template FP_NR<long double> svp_probability<FP_NR<long double>>(const vector<double> &pr);

#endif

#ifdef FPLLL_WITH_QD

template class Pruner<FP_NR<dd_real>>;
template void prune<FP_NR<dd_real>>(Pruning &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<dd_real>>(Pruning &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<dd_real> svp_probability<FP_NR<dd_real>>(const Pruning &pruning);
template FP_NR<dd_real> svp_probability<FP_NR<dd_real>>(const vector<double> &pr);

template class Pruner<FP_NR<qd_real>>;
template void prune<FP_NR<qd_real>>(Pruning &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<qd_real>>(Pruning &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<qd_real> svp_probability<FP_NR<qd_real>>(const Pruning &pruning);
template FP_NR<qd_real> svp_probability<FP_NR<qd_real>>(const vector<double> &pr);

#endif

#ifdef FPLLL_WITH_DPE

template class Pruner<FP_NR<dpe_t>>;
template void prune<FP_NR<dpe_t>>(Pruning &,const double, const double, const vector<double> &, const double, const PrunerMetric, const int);
template void prune<FP_NR<dpe_t>>(Pruning &,const double, const double, const vector<vector<double>> &, const double, const PrunerMetric, const int);
template FP_NR<dpe_t> svp_probability<FP_NR<dpe_t>>(const Pruning &pruning);
template FP_NR<dpe_t> svp_probability<FP_NR<dpe_t>>(const vector<double> &pr);

#endif
/* clang-format on */

FPLLL_END_NAMESPACE
