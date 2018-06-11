#include "fplll.h"

FPLLL_BEGIN_NAMESPACE

#define OPTIMIZE_PROB_MINSTEP 1e-4
#define OPTIMIZE_PROB_MAXSTEP 1e4

template <class FT> void Pruner<FT>::optimize_coefficients_prob_incr(/*io*/ vector<double> &pr)
{
  int dn = pr.size();
  int tours;
  double normalized;
  FT old_c0, old_c1, old_prob, old_cfs;
  vec b(dn), old_b(dn), old_b2(dn);
  vector<double> detailed_cost(dn);
  vector<double> weight(dn);
  bool not_changed;

  load_coefficients(b, pr);

  // incr b until achieve target or quit if too much or quite if impossible
  tours = 0;
  while (1)
  {
    if (tours > OPTIMIZE_PROB_MAXSTEP)
      break;
    tours++;

    old_prob = measure_metric(b);
    if (old_prob >= target)
      break;

    old_cfs    = single_enum_cost(b, &(detailed_cost));
    normalized = 0.0;
    for (int i = 0; i < dn; i++)
    {
      weight[i] = 0.0;
      for (int j = i; j < dn; j++)
      {
        weight[i] = weight[i] + detailed_cost[j];
      }
      weight[i] = 1.0 / weight[i];
      if (weight[i] < OPTIMIZE_PROB_MINSTEP)
        weight[i] = OPTIMIZE_PROB_MINSTEP;
      normalized += weight[i];
    }
    for (int i = 0; i < dn; i++)
    {
      weight[i] = weight[i] / normalized;
    }
    for (int i = dn - 1; i >= 0; --i)
    {
      old_b[i] = b[i];
      b[i]     = b[i] + weight[i];
      if (b[i] >= 1.0)
        b[i] = 1.0;
    }

    enforce(b);

    not_changed = true;
    for (int i = dn - 1; i >= 0; --i)
    {
      if (b[i] != old_b[i])
        not_changed = false;
    }
    if (not_changed)
      break;
  }

  save_coefficients(pr, b);
}

template <class FT> void Pruner<FT>::optimize_coefficients_prob_decr(/*io*/ vector<double> &pr)
{
  int dn = pr.size();
  int tours;
  double normalized;
  FT old_c0, old_c1, old_prob, old_cfs;
  vec b(dn), old_b(dn), old_b2(dn);
  vector<double> detailed_cost(dn);
  vector<double> weight(dn);
  bool not_changed;

  load_coefficients(b, pr);
  // decr b until achieve target
  tours = 0;
  while (1)
  {
    if (tours > OPTIMIZE_PROB_MAXSTEP)
      break;
    tours++;

    old_prob = measure_metric(b);
    if (old_prob <= target)
      break;

    old_cfs    = single_enum_cost(b, &(detailed_cost));
    normalized = 0.0;
    for (int i = 0; i < dn; i++)
    {
      weight[i] = 0.0;
      for (int j = i; j < dn; j++)
      {
        weight[i] = weight[i] + detailed_cost[j];
      }
      weight[i] = 1.0 / weight[i];
      if (weight[i] < OPTIMIZE_PROB_MINSTEP)
        weight[i] = OPTIMIZE_PROB_MINSTEP;
      normalized += weight[i];
    }
    for (int i = 0; i < dn; i++)
    {
      weight[i] = weight[i] / normalized;
      // cout << weight[i] << " ";
    }

    for (int i = dn - 1; i >= 0; --i)
    {
      old_b[i] = b[i];
      b[i]     = b[i] - weight[i];
      if (b[i] < OPTIMIZE_PROB_MINSTEP)
        b[i] = OPTIMIZE_PROB_MINSTEP;
    }

    enforce(b);

    not_changed = true;
    for (int i = dn - 1; i >= 0; --i)
    {
      if (b[i] != old_b[i])
        not_changed = false;
    }
    if (not_changed)
    {
      break;
    }
  }
  save_coefficients(pr, b);
}

template <class FT> void Pruner<FT>::optimize_coefficients_prob_tune(/*io*/ vector<double> &pr)
{
  int dn = pr.size();
  int tours;
  FT prob, ratio;
  vec b(dn), old_b(dn), old_b2(dn);
  vector<double> detailed_cost(dn);
  vector<double> weight(dn);
  bool not_changed;

  load_coefficients(b, pr);

  // incr b until achieve target
  tours = 0;
  while (1)
  {
    tours++;

    prob  = measure_metric(b);
    ratio = prob / target;

    // good enough
    if (ratio < 1.05 && ratio > 0.95)
      break;

    // tune
    if (ratio < 1)
    {
      for (int i = dn - 1; i >= 0; --i)
      {
        old_b[i] = b[i];
        b[i]     = b[i] + OPTIMIZE_PROB_MINSTEP;
        if (b[i] >= 1.0)
          b[i] = 1.0;
      }
    }
    else
    {
      for (int i = dn - 1; i >= 0; --i)
      {
        old_b[i] = b[i];
        b[i]     = b[i] - OPTIMIZE_PROB_MINSTEP;
        if (b[i] < OPTIMIZE_PROB_MINSTEP)
          b[i] = OPTIMIZE_PROB_MINSTEP;
      }
    }

    enforce(b);

    not_changed = true;
    for (int i = dn - 1; i >= 0; --i)
    {
      if (b[i] != old_b[i])
        not_changed = false;
    }
    if (not_changed)
      break;
  }

  save_coefficients(pr, b);
}

FPLLL_END_NAMESPACE
