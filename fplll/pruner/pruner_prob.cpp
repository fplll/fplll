#include "fplll.h"

FPLLL_BEGIN_NAMESPACE

template <class FT> inline FT Pruner<FT>::svp_probability_evec(/*i*/ const evec &b)
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
  FT res  = dvol / (dxn - 1.);

  if (!res.is_finite())
  {
    throw std::range_error("NaN or inf in svp_probability");
  }
  return res;
}

/**
 * hopefully lower bound for bounding b
 */
template <class FT> FT Pruner<FT>::svp_probability_lower(/*i*/ const vec &b)
{
  evec b_lower(d);
  for (int i = 0; i < d; ++i)
  {
    b_lower[i] = b[2 * i];
  }
  return svp_probability_evec(b_lower);
}

/**
 * upper bound for bounding b
 */
template <class FT> FT Pruner<FT>::svp_probability_upper(/*i*/ const vec &b)
{
  evec b_upper(d);
  for (int i = 0; i < d; ++i)
  {
    b_upper[i] = b[2 * i + 1];
  }
  return svp_probability_evec(b_upper);
}

/**
 * upper estimate for bounding b
 */
template <class FT> inline FT Pruner<FT>::svp_probability(/*i*/ const vec &b)
{
  /*
  if (!shape_loaded)
  {
    throw std::invalid_argument("No basis shape was loaded");
  }
  */
  if (b.size() == (unsigned int)d)
  {
    return svp_probability_evec(b);
  }
  else
  {
    FT pl = svp_probability_lower(b);
    FT pu = svp_probability_upper(b);
    return (pl + pu) / 2.0;
  }
}

template <class FT> inline FT Pruner<FT>::expected_solutions_evec(/*i*/ const evec &b)
{
  int j  = d * 2 - 1;
  FT tmp = log(relative_volume(d, b));
  tmp += log(tabulated_ball_vol[j + 1]);
  tmp += (log(normalized_radius) + log(b[j / 2]) / 2.0) * (j + 1);
  tmp += log(ipv[j]);
  tmp += log(symmetry_factor);
  tmp = exp(tmp);

  if (!tmp.is_finite())
  {
    throw std::range_error("Error: NaN or inf in expected_solutions");
  }
  return tmp;
}

template <class FT> FT Pruner<FT>::expected_solutions_lower(/*i*/ const vec &b)
{
  evec b_lower(d);
  for (int i = 0; i < d; ++i)
  {
    b_lower[i] = b[2 * i];
  }
  return expected_solutions_evec(b_lower);
}

template <class FT> FT Pruner<FT>::expected_solutions_upper(/*i*/ const vec &b)
{
  evec b_upper(d);
  for (int i = 0; i < d; ++i)
  {
    b_upper[i] = b[2 * i + 1];
  }
  return expected_solutions_evec(b_upper);
}

template <class FT> inline FT Pruner<FT>::expected_solutions(/*i*/ const vec &b)
{
  if (!shape_loaded)
  {
    throw std::invalid_argument("No basis shape was loaded");
  }
  if (b.size() == (unsigned int)d)
  {
    return expected_solutions_evec(b);
  }
  else
  {
    FT pl = expected_solutions_lower(b);
    FT pu = expected_solutions_upper(b);
    return (pl + pu) / 2.0;
  }
}

template <class FT> inline FT Pruner<FT>::measure_metric(/*i*/ const vec &b)
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

FPLLL_END_NAMESPACE
