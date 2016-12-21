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

template <class FT> FT svp_probability(const Pruning &pruning)
{
  Pruner<FT> pru;
  return pru.measure_metric(pruning.coefficients);
}

template <class FT> FT svp_probability(const vector<double> &pr)
{
  Pruner<FT> pru;
  return pru.measure_metric(pr);
}

template <class FT> void Pruner<FT>::set_tabulated_consts()
{
  for (int i = 0; i < PRUNER_MAX_N; ++i)
  {
    tabulated_factorial[i] = pre_factorial[i];
    tabulated_ball_vol[i]  = pre_ball_vol[i];
  }
  return;
}

/// PUBLIC METHODS
template <class FT>
void Pruner<FT>::load_basis_shape(const vector<double> &gso_sq_norms, bool reset_renorm)
{
  n = gso_sq_norms.size();
  d = n / 2;
  if (!d)
  {
    throw std::runtime_error("Inside Pruner : Needs a dimension n>1");
  }
  FT logvol, tmp;
  logvol = 0.0;
  for (size_t i = 0; i < n; ++i)
  {
    r[i] = gso_sq_norms[n - 1 - i];
    logvol += log(r[i]);
  }
  if (reset_renorm)
  {
    renormalization_factor = exp(logvol / (-1.0 * n));
  }

  for (size_t i = 0; i < n; ++i)
  {
    r[i] *= renormalization_factor;
  }
  tmp = 1.;
  for (size_t i = 0; i < 2 * d; ++i)
  {
    tmp *= sqrt(r[i]);
    ipv[i] = 1.0 / tmp;
  }
}

template <class FT>
void Pruner<FT>::load_basis_shapes(const vector<vector<double>> &gso_sq_norms_vec)
{
  vec sum_ipv;
  n = gso_sq_norms_vec[0].size();
  for (size_t i = 0; i < n; ++i)
  {
    sum_ipv[i] = 0.;
  }
  int count = gso_sq_norms_vec.size();
  for (int k = 0; k < count; ++k)
  {
    if (gso_sq_norms_vec[k].size() != n)
    {
      throw std::runtime_error("Inside Pruner : loading several bases with different dimensions");
    }
    bool reset_renorm = (k == 0);
    load_basis_shape(gso_sq_norms_vec[k], reset_renorm);
    for (size_t i = 0; i < n; ++i)
    {
      sum_ipv[i] += ipv[i];
    }
  }
  for (size_t i = 0; i < n; ++i)
  {
    ipv[i] = sum_ipv[i] / (1.0 * count);
  }
}
template <class FT>
void Pruner<FT>::optimize_coefficients(/*io*/ vector<double> &pr, /*i*/ const bool reset)
{
  evec b;
  for (int i = 0; i < PRUNER_MAX_D; ++i)
  {
    b[i] = 0.;
  }
  init_coefficients(b);
  if (!reset)
  {
    load_coefficients(b, pr);
  }
  descent(b);
  save_coefficients(pr, b);
}

// PRIVATE METHODS

template <class FT>
void Pruner<FT>::load_coefficients(/*o*/ evec &b, /*i*/ const vector<double> &pr)
{
  for (size_t i = 0; i < d; ++i)
  {
    b[i] = pr[n - 1 - 2 * i];
  }
  if (enforce_bounds(b))
  {
    throw std::runtime_error("Inside Pruner : Ill formed pruning coefficients (must be decreasing, "
                             "starting with two 1.0)");
  }
}

template <class FT> int Pruner<FT>::check_basis_loaded()
{
  if (d)
  {
    return 0;
  }
  throw std::runtime_error("Inside Pruner : No basis loaded");
  return 1;
}

template <class FT>
void Pruner<FT>::save_coefficients(/*o*/ vector<double> &pr, /*i*/ const evec &b)
{
  pr.resize(n);
  for (size_t i = 0; i < d; ++i)
  {
    pr[n - 1 - 2 * i] = b[i].get_d();
    pr[n - 2 - 2 * i] = b[i].get_d();
  }
  pr[0] = 1.;
}

template <class FT> inline int Pruner<FT>::enforce_bounds(/*io*/ evec &b, /*opt i*/ const int j)
{
  int status = 0;
  if ((b[d - 1] < .999) & (d - j != 1))
  {
    status   = 1;
    b[d - 1] = 1.;
  }
  for (size_t i = 0; i < d; ++i)
  {
    if (b[i] > 1.0001)
    {
      status = 1;
    }
    if (b[i] > 1)
    {
      b[i] = 1.0;
    }
    if (b[i] <= .05)
      b[i] = .05;
  }
  for (size_t i = j; i < d - 1; ++i)
  {
    if (b[i + 1] < b[i])
    {
      if (b[i + 1] + .001 < b[i])
      {
        status = 1;
      }
      b[i + 1] = b[i];
    }
  }
  for (int i = j - 1; i >= 0; --i)
  {
    if (b[i + 1] < b[i])
    {
      if (b[i + 1] + .001 < b[i])
      {
        status = 1;
      }
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
  poly P;
  P[0]   = 1;
  int ld = 0;
  for (int i = rd - 1; i >= 0; --i)
  {
    integrate_poly(ld, P);
    ld++;
    P[0] = -1.0 * eval_poly(ld, P, b[i] / b[rd - 1]);
  }
  if (rd % 2)
  {
    return -1.0 * P[0] * tabulated_factorial[rd];
  }
  else
  {
    return P[0] * tabulated_factorial[rd];
  }
}

template <class FT>
inline FT Pruner<FT>::single_enum_cost(/*i*/ const evec &b, vector<double> *detailed_cost)
{

  if (detailed_cost)
  {
    detailed_cost->resize(n);
  }
  vec rv;  // Relative volumes at each level

  for (size_t i = 0; i < d; ++i)
  {

    rv[2 * i + 1] = relative_volume(i + 1, b);
  }

  rv[0] = 1;
  for (size_t i = 1; i < d; ++i)
  {
    rv[2 * i] = sqrt(rv[2 * i - 1] * rv[2 * i + 1]);  // Interpolate even values
  }

  FT total;
  total = 0.0;
  FT normalized_radius;
  normalized_radius = sqrt(enumeration_radius * renormalization_factor);

  FT normalized_radius_pow = normalized_radius;
  for (size_t i = 0; i < 2 * d; ++i)
  {
    FT tmp;

    tmp = normalized_radius_pow * rv[i] * tabulated_ball_vol[i + 1] *
          sqrt(pow_si(b[i / 2], 1 + i)) * ipv[i];

    if (detailed_cost)
    {
      (*detailed_cost)[i] = tmp.get_d();
    }

    total += tmp;
    normalized_radius_pow *= normalized_radius;
  }
  total /= symmetry_factor;
  return total;
}

template <class FT> inline FT Pruner<FT>::svp_probability(/*i*/ const evec &b)
{
  evec b_minus_db;
  FT dx = shell_ratio;
  for (size_t i = 0; i < d; ++i)
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
  FT normalized_radius;
  normalized_radius = sqrt(enumeration_radius * renormalization_factor);

  FT vol = relative_volume(d, b);
  vol *= tabulated_ball_vol[2 * d - 1];
  vol *= pow_si(normalized_radius, 2 * d);
  vol *= ipv[2 * d - 1];

  return vol;
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
  evec bpDb;
  res[d - 1] = 0.0;
  for (size_t i = 0; i < d - 1; ++i)
  {
    bpDb = b;
    bpDb[i] *= (1.0 - epsilon);
    enforce_bounds(bpDb, i);
    FT X = repeated_enum_cost(bpDb);

    bpDb = b;
    bpDb[i] *= (1.0 + epsilon);
    enforce_bounds(bpDb, i);
    FT Y   = repeated_enum_cost(bpDb);
    res[i] = (log(X) - log(Y)) / epsilon;
  }
}

template <class FT> int Pruner<FT>::improve(/*io*/ evec &b)
{

  FT cf     = repeated_enum_cost(b);
  FT old_cf = cf;
  evec newb;
  evec gradient;
  repeated_enum_cost_gradient(b, gradient);
  FT norm = 0.0;

  // normalize the gradient
  for (size_t i = 0; i < d; ++i)
  {
    norm += gradient[i] * gradient[i];
    newb[i] = b[i];
  }

  norm /= (1.0 * (1. * d));
  norm = sqrt(norm);
  if (norm <= 0.)
    return 0;

  for (size_t i = 0; i < d; ++i)
  {
    gradient[i] /= norm;
  }
  FT new_cf;

  FT step = min_step;
  size_t i;

  for (i = 0;; ++i)
  {
    for (size_t i = 0; i < d; ++i)
    {
      newb[i] = newb[i] + step * gradient[i];
    }

    enforce_bounds(newb);
    new_cf = repeated_enum_cost(newb);

    if (new_cf >= cf)
    {
      break;
    }
    b  = newb;
    cf = new_cf;
    step *= step_factor;
  }
  if (cf > old_cf * min_cf_decrease)
  {
    return 0;
  }
  return i;
}

template <class FT> void Pruner<FT>::init_coefficients(evec &b)
{
  for (size_t i = 0; i < d; ++i)
  {
    b[i] = .1 + ((1. * i) / d);
  }
  enforce_bounds(b);
}

template <class FT> void Pruner<FT>::descent(/*io*/ evec &b)
{
  if (method == PRUNER_METHOD_GREEDY)
  {
    greedy(b);
    return;
  }

  if ((method == PRUNER_METHOD_GRADIENT) || (method == PRUNER_METHOD_HYBRID))
  {
    while (improve(b))
    {
    };
  };
  if ((method == PRUNER_METHOD_NM) || (method == PRUNER_METHOD_HYBRID))
  {
    while (nelder_mead(b))
    {
    };
  };
}

template <class FT> void Pruner<FT>::greedy(evec &b)
{
  for (size_t i = 0; i < d; ++i)
  {
    b[i] = 1.;
  }
  enforce_bounds(b);

  evec newb;

  FT normalized_radius;
  normalized_radius = sqrt(enumeration_radius * renormalization_factor);

  FT min, max, val, tmp, goal;
  if (verbosity)
  {
    cerr << "Starting Greedy pruning" << endl;
  }
  for (size_t j = 1; j < 2 * d; j += 2)
  {
    val = 1.;
    max = 1.;
    min = 0.025;
    if (j == 2 * d - 1)
    {
      goal = target;
    }
    else
    {
      goal = preproc_cost / (2 * d);
    }
    int count = 0;
    tmp       = 0.;
    while ((count < 12) && (min < .99))
    {
      if (val < .05)
      {
        enumeration_radius /= 2.;
        greedy(b);
        return;
      }
      count++;
      newb        = b;
      newb[j / 2] = val;
      enforce_bounds(newb, j / 2);

      tmp = relative_volume((j + 1) / 2, newb);
      tmp *= tabulated_ball_vol[j + 1];
      tmp *= pow_si(normalized_radius * sqrt(newb[j / 2]), j + 1);
      tmp *= ipv[j];

      if (tmp > goal)
      {
        max = val;
      }
      else
      {
        min = val;
      }
      val = (min + max) / 2.;
    }
    if (verbosity)
    {
      cerr << j << " : " << val << " ~ " << tmp.get_d() << " G " << goal << endl;
    }
    b[j / 2] = val;
    enforce_bounds(b, j / 2);
  }

  FT factor = b[d - 1];
  for (size_t i = 0; i < d; ++i)
  {
    b[i] /= factor;
  }
  enforce_bounds(b);
  enumeration_radius *= factor;
  normalized_radius = sqrt(enumeration_radius * renormalization_factor);

  tmp = relative_volume(d, b);
  tmp *= tabulated_ball_vol[2 * d - 1];
  tmp *= pow_si(normalized_radius * sqrt(b[d - 1]), 2 * d);
  tmp *= ipv[2 * d - 1];
}

// Nelder-Mead method. Following the notation of
// https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method

#define ND_ALPHA 1
#define ND_GAMMA 2
#define ND_RHO 0.5
#define ND_SIGMA 0.5
#define ND_INIT_WIDTH 0.05

template <class FT> int Pruner<FT>::nelder_mead(/*io*/ evec &b)
{
  size_t l = d + 1;
  evec *bs = new evec[l];  // The simplexe (d+1) vector of dim d
  FT *fs   = new FT[l];    // Values of f at the simplex vertices

  for (size_t i = 0; i < l; ++i)  // Intialize the simplex
  {
    bs[i] = b;  // One of the is b
    if (i < d)
    {
      if (bs[i][i] < .5)
      {
        bs[i][i] += ND_INIT_WIDTH;  // the other are perturbations on b of +/- const
      }
      else
      {
        bs[i][i] -=
            ND_INIT_WIDTH;  // sign is chosen to avoid getting close to the border of the domain
      }
    }
    enforce_bounds(bs[i]);
    fs[i] = repeated_enum_cost(bs[i]);  // initialize the value
  }

  FT init_cf = fs[l - 1];

  evec bo;  // centeroid
  FT fo;    // value at the centroid

  FT fs_maxi_last;  // value of the last centroid

  if (verbosity)
  {
    cerr << "  Starting nelder_mead cf = " << init_cf << " proba = " << svp_probability(b) << endl;
  }
  unsigned int counter = 0;
  size_t mini = 0, maxi = 0, maxi2 = 0;
  while (1)  // Main loop
  {
    mini = maxi = maxi2 = 0;
    for (size_t i = 0; i < d; ++i)
      bo[i]       = bs[0][i];
    ////////////////
    // step 1. and 2. : Order and centroid
    ////////////////
    for (size_t i = 1; i < l; ++i)  // determine min and max, and centroid
    {
      if (fs[i] < fs[mini])
        mini = i;
      if (fs[i] > fs[maxi])
        maxi = i;
      for (size_t j = 0; j < d; ++j)
        bo[j] += bs[i][j];
    }
    FT tmp;
    tmp = l;
    for (size_t i = 0; i < d; ++i)
      bo[i] /= tmp;  // Centroid calculated

    if (!counter)
      fs_maxi_last = fs[maxi];

    if (!maxi)
      maxi2++;
    for (size_t i = 1; i < l; ++i)  // determine min and max, and centroid
    {
      if ((fs[i] > fs[maxi2]) && (i != maxi))
        maxi2 = i;
    }

    if (enforce_bounds(bo))
    {
      throw std::runtime_error("Concavity says that should not happen.");
    }

    if (verbosity)
    {
      cerr << "  melder_mead step " << counter << "cf = " << fs[mini]
           << " proba = " << measure_metric(bs[mini]) << " cost = " << single_enum_cost(bs[mini])
           << endl;
      for (size_t i = 0; i < d; ++i)
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

    for (size_t i = 0; i < l; ++i)  // determine second best
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

    evec br;  // reflected point
    FT fr;    // Value at the reflexion point
    for (size_t i = 0; i < d; ++i)
      br[i]       = bo[i] + ND_ALPHA * (bo[i] - bs[maxi][i]);
    enforce_bounds(br);
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
      evec be;
      FT fe;
      for (size_t i = 0; i < d; ++i)
        be[i]       = bo[i] + ND_GAMMA * (br[i] - bo[i]);
      enforce_bounds(be);
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

    evec bc;
    FT fc;
    for (size_t i = 0; i < d; ++i)
      bc[i]       = bo[i] + ND_RHO * (bs[maxi][i] - bo[i]);
    enforce_bounds(bc);
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
    for (size_t j = 0; j < l; ++j)
    {
      for (size_t i = 0; i < d; ++i)
      {
        bs[j][i] = bs[mini][i] + ND_SIGMA * (bs[j][i] - bs[mini][i]);
      }
      enforce_bounds(bs[j]);
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

  delete[] fs;
  delete[] bs;

  return improved;  // Has MN made any progress
}

template <class FT>
void prune(/*output*/ Pruning &pruning,
           /*inputs*/ double &enumeration_radius, const double preproc_cost, const double target,
           vector<double> &r, const PrunerMethod method, const PrunerMetric metric, bool reset)
{
  Pruner<FT> pruner(enumeration_radius, preproc_cost, target, method, metric);
  pruner.load_basis_shape(r);
  pruner.optimize_coefficients(pruning.coefficients);
  pruner.single_enum_cost(pruning.coefficients, &(pruning.detailed_cost));
  enumeration_radius  = pruner.enumeration_radius.get_d();
  pruning.metric      = metric;
  pruning.expectation = pruner.measure_metric(pruning.coefficients);
}

template <class FT>
void prune(/*output*/ Pruning &pruning,
           /*inputs*/ double &enumeration_radius, const double preproc_cost, const double target,
           vector<vector<double>> &rs, const PrunerMethod method, const PrunerMetric metric,
           bool reset)
{
  Pruner<FT> pruner(enumeration_radius, preproc_cost, target, method, metric);
  pruner.load_basis_shapes(rs);
  pruner.optimize_coefficients(pruning.coefficients);
  pruner.single_enum_cost(pruning.coefficients, &(pruning.detailed_cost));
  enumeration_radius  = pruner.enumeration_radius.get_d();
  pruning.metric      = metric;
  pruning.expectation = pruner.measure_metric(pruning.coefficients);
}

/** instantiate functions **/

//////////////////////////////////////// Instantiations

// DOUBLE

template class Pruner<FP_NR<double>>;

template void prune<FP_NR<double>>(Pruning &, double &, const double, const double,
                                   vector<double> &, const PrunerMethod, const PrunerMetric, bool);
template void prune<FP_NR<double>>(Pruning &, double &, const double, const double,
                                   vector<vector<double>> &, const PrunerMethod, const PrunerMetric,
                                   bool);
template FP_NR<double> svp_probability<FP_NR<double>>(const Pruning &pruning);
template FP_NR<double> svp_probability<FP_NR<double>>(const vector<double> &pr);

template class Pruner<FP_NR<mpfr_t>>;

template void prune<FP_NR<mpfr_t>>(Pruning &, double &, const double, const double,
                                   vector<double> &, const PrunerMethod, const PrunerMetric, bool);
template void prune<FP_NR<mpfr_t>>(Pruning &, double &, const double, const double,
                                   vector<vector<double>> &, const PrunerMethod, const PrunerMetric,
                                   bool);
template FP_NR<mpfr_t> svp_probability<FP_NR<mpfr_t>>(const Pruning &pruning);
template FP_NR<mpfr_t> svp_probability<FP_NR<mpfr_t>>(const vector<double> &pr);

// LD

#ifdef FPLLL_WITH_LONG_DOUBLE

template class Pruner<FP_NR<long double>>;
template void prune<FP_NR<long double>>(Pruning &, double &, const double, const double,
                                        vector<double> &, const PrunerMethod, const PrunerMetric,
                                        bool);
template void prune<FP_NR<long double>>(Pruning &, double &, const double, const double,
                                        vector<vector<double>> &, const PrunerMethod,
                                        const PrunerMetric, bool);
template FP_NR<long double> svp_probability<FP_NR<long double>>(const Pruning &pruning);
template FP_NR<long double> svp_probability<FP_NR<long double>>(const vector<double> &pr);

#endif

#ifdef FPLLL_WITH_QD

template class Pruner<FP_NR<dd_real>>;
template void prune<FP_NR<dd_real>>(Pruning &, double &, const double, const double,
                                    vector<double> &, const PrunerMethod, const PrunerMetric, bool);
template void prune<FP_NR<dd_real>>(Pruning &, double &, const double, const double,
                                    vector<vector<double>> &, const PrunerMethod,
                                    const PrunerMetric, bool);
template FP_NR<dd_real> svp_probability<FP_NR<dd_real>>(const Pruning &pruning);
template FP_NR<dd_real> svp_probability<FP_NR<dd_real>>(const vector<double> &pr);

template class Pruner<FP_NR<qd_real>>;
template void prune<FP_NR<qd_real>>(Pruning &, double &, const double, const double,
                                    vector<double> &, const PrunerMethod, const PrunerMetric, bool);
template void prune<FP_NR<qd_real>>(Pruning &, double &, const double, const double,
                                    vector<vector<double>> &, const PrunerMethod,
                                    const PrunerMetric, bool);
template FP_NR<qd_real> svp_probability<FP_NR<qd_real>>(const Pruning &pruning);
template FP_NR<qd_real> svp_probability<FP_NR<qd_real>>(const vector<double> &pr);

#endif

#ifdef FPLLL_WITH_DPE

template class Pruner<FP_NR<dpe_t>>;
template void prune<FP_NR<dpe_t>>(Pruning &, double &, const double, const double, vector<double> &,
                                  const PrunerMethod, const PrunerMetric, bool);
template void prune<FP_NR<dpe_t>>(Pruning &, double &, const double, const double,
                                  vector<vector<double>> &, const PrunerMethod, const PrunerMetric,
                                  bool);
template FP_NR<dpe_t> svp_probability<FP_NR<dpe_t>>(const Pruning &pruning);
template FP_NR<dpe_t> svp_probability<FP_NR<dpe_t>>(const vector<double> &pr);

#endif

FPLLL_END_NAMESPACE
