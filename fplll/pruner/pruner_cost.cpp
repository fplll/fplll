#include "fplll.h"

FPLLL_BEGIN_NAMESPACE

/**
 * single enumation cost using half-vector coefficients b
 */
template <class FT>
inline FT Pruner<FT>::single_enum_cost_evec(/*i*/ const evec &b,
                                            /*o*/ vector<double> *detailed_cost)
{
  if (!shape_loaded)
  {
    throw std::invalid_argument("Error: No basis shape was loaded");
  }

  if (detailed_cost)
  {
    detailed_cost->resize(n);
  }
  vec rv(n);  // Relative volumes at each level

  /* even dimension of 2*i, e.g. 2, 4, */
  for (int i = 0; i < d; ++i)
  {
    rv[2 * i + 1] = relative_volume(i + 1, b);
  }

  /* old dimension, e.g. get dim 3 by project (dim 2 and 4) */
  rv[0] = 1;
  for (int i = 1; i < d; ++i)
  {
    rv[2 * i] = sqrt(rv[2 * i - 1] * rv[2 * i + 1]);  // Interpolate
                                                      // even values
    // another way
    // rv[2 * i] = (rv[2 * i - 1] + rv[2 * i + 1]) / 2.0;
  }

  FT total;
  total                    = 0.0;
  FT normalized_radius_pow = normalized_radius;

  for (int i = 0; i < 2 * d; ++i)
  {
    FT tmp;

    // ipv is normalized and larger
    // normalized_radius_pow = (R/det^1/n)^i and normalized to be
    // smaller, where R is the actual searching radius in SVP..
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
  if (!total.is_finite())
  {
    throw std::range_error("NaN or inf in single_enum_cost");
  }
  return total;
}

/**
 * lower bound of single enumation cost
 */
template <class FT>
inline FT Pruner<FT>::single_enum_cost_lower(/*i*/ const vec &b, vector<double> *detailed_cost)
{
  evec b_lower(d);
  for (int i = 0; i < d; ++i)
  {
    b_lower[i] = b[2 * i];
  }
  return single_enum_cost_evec(b_lower, detailed_cost);
}

/**
 * upper bound of single enumation cost
 */
template <class FT>
inline FT Pruner<FT>::single_enum_cost_upper(/*i*/ const vec &b, vector<double> *detailed_cost)
{
  evec b_upper(d);
  for (int i = 0; i < d; ++i)
  {
    b_upper[i] = b[2 * i + 1];
  }
  return single_enum_cost_evec(b_upper, detailed_cost);
}

/**
 * single enumation cost
 */
template <class FT>
inline FT Pruner<FT>::single_enum_cost(/*i*/ const vec &b, vector<double> *detailed_cost)
{
  if (b.size() == (unsigned int)d)
  {
    return single_enum_cost_evec(b, detailed_cost);
  }
  else
  {
    FT cl = single_enum_cost_lower(b, detailed_cost);
    FT cu = single_enum_cost_upper(b, detailed_cost);
    return (cl + cu) / 2.0;
  }
}

template <class FT> void Pruner<FT>::target_function_gradient(/*i*/ const vec &b, /*o*/ vec &res)
{

  int dn = b.size();
  vec b_plus_db(dn);
  res[dn - 1] = 0.0;  // Force null gradient on the last coordinate : don't touch this coeff
  for (int i = 0; i < dn - 1; ++i)
  {
    b_plus_db = b;
    b_plus_db[i] *= (1.0 - epsilon);
    enforce(b_plus_db, i);
    FT X = target_function(b_plus_db);

    b_plus_db = b;
    b_plus_db[i] *= (1.0 + epsilon);
    enforce(b_plus_db, i);
    FT Y   = target_function(b_plus_db);
    res[i] = (log(X) - log(Y)) / epsilon;
  }
}

template <class FT> inline FT Pruner<FT>::target_function(/*i*/ const vec &b)
{
  if (metric == PRUNER_METRIC_PROBABILITY_OF_SHORTEST)
  {
    FT probability = svp_probability(b);
    FT trials      = log(1.0 - target) / log(1.0 - probability);
    if (!trials.is_finite())
    {
      throw std::range_error("NaN or inf in target_function (METRIC_PROBABILITY_OF_SHORTEST)");
    }
    trials = trials < 1.0 ? 1.0 : trials;
    return single_enum_cost(b) * trials + preproc_cost * (trials - 1.0);
  }
  else if (metric == PRUNER_METRIC_EXPECTED_SOLUTIONS)
  {
    FT expected = expected_solutions(b);
    FT trials   = target / expected;
    if (!trials.is_finite())
    {
      throw std::range_error("NaN or inf in target_function (METRIC_EXPECTED_SOLUTION)");
    }
    // if expected solutions > 1, set trial = 1
    trials = trials < 1.0 ? 1.0 : trials;
    return single_enum_cost(b) * trials + preproc_cost * (trials - 1.0);
  }
  else
  {
    throw std::invalid_argument("Pruner was set to an unknown metric");
  }
}

template <class FT> inline FT Pruner<FT>::repeated_enum_cost(/*i*/ const vec &b)
{
  if (metric == PRUNER_METRIC_PROBABILITY_OF_SHORTEST)
  {
    FT probability = svp_probability(b);
    FT trials      = log(1.0 - target) / log(1.0 - probability);
    if (!trials.is_finite())
    {
      throw std::range_error("NaN or inf in repeated_enum_cost (METRIC_PROBABILITY_OF_SHORTEST)");
    }
    trials = trials < 1.0 ? 1.0 : trials;
    return single_enum_cost(b) * trials + preproc_cost * (trials - 1.0);
  }
  else if (metric == PRUNER_METRIC_EXPECTED_SOLUTIONS)
  {
    FT expected = expected_solutions(b);
    FT trials   = 1.0 / expected;
    if (!trials.is_finite())
    {
      throw std::range_error("NaN or inf in repeated_enum_cost (METRIC_EXPECTED_SOLUTION)");
    }
    // if expected solutions > 1, set trial = 1
    trials = trials < 1.0 ? 1.0 : trials;
    return single_enum_cost(b) * trials + preproc_cost * (trials - 1.0);
  }
  else
  {
    throw std::invalid_argument("Pruner was set to an unknown metric");
  }
}

FPLLL_END_NAMESPACE
