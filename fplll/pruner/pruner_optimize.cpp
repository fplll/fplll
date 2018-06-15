#include "fplll.h"

FPLLL_BEGIN_NAMESPACE

//#define DEBUG_PRUNER_OPTIMIZE
#define NUM_OPTIMIZATION_TOURS 3

template <class FT> void Pruner<FT>::optimize_coefficients_cost(/*io*/ vector<double> &pr)
{

  FT old_c0, old_c1, new_c, min_c;
  vec b(n), best_b(n);

  // step 1 use half coefficients only
  optimize_coefficients_evec(pr);

  load_coefficients(b, pr);
  old_c0 = repeated_enum_cost(b);
  min_c  = old_c0;

#ifdef DEBUG_PRUNER_OPTIMIZE
  FT old_sc   = single_enum_cost(b);
  FT old_prob = measure_metric(b);
  cerr << "# [Stage ho] half-optimization done." << endl;
  cerr << "# [Stage ho]    all_enum_cost = " << old_c0 << endl;
  cerr << "# [Stage ho] single_enum_cost = " << old_sc << endl;
  cerr << "# [Stage ho] succ_probability = " << old_prob << endl;
#endif

  // step 2 use full coefficients if enabled
  if (flags & PRUNER_FULL)
  {
    int tours = 0;
    while (1)
    {
      tours++;

      // initial cost
      load_coefficients(b, pr);
      old_c0 = repeated_enum_cost(b);

      // step 2.1 full tuning
      optimize_coefficients_tune_single(pr);
      optimize_coefficients_tune_prob(pr);
      optimize_coefficients_smooth(pr);

      // update best after tuning
      load_coefficients(b, pr);
      old_c1 = repeated_enum_cost(b);
      if (old_c1 < min_c)
      {
        min_c  = old_c1;
        best_b = b;
      }

#ifdef DEBUG_PRUNER_OPTIMIZE
      old_sc   = single_enum_cost(b);
      old_prob = measure_metric(b);
      cerr << "# [Stage ft] full-tuning done." << endl;
      cerr << "# [Stage ft]    all_enum_cost = " << old_c1 << endl;
      cerr << "# [Stage ft] single_enum_cost = " << old_sc << endl;
      cerr << "# [Stage ft] succ_probability = " << old_prob << endl;
#endif

      // step 2.2 full optimization
      optimize_coefficients_full(pr);

      load_coefficients(b, pr);
      new_c = repeated_enum_cost(b);
      if (new_c < min_c)
      {
        min_c  = new_c;
        best_b = b;
      }

#ifdef DEBUG_PRUNER_OPTIMIZE
      old_sc   = single_enum_cost(b);
      old_prob = measure_metric(b);
      cerr << "# [Stage fo] full-optimization done." << endl;
      cerr << "# [Stage fo]    all_enum_cost = " << new_c << endl;
      cerr << "# [Stage fo] single_enum_cost = " << old_sc << endl;
      cerr << "# [Stage fo] succ_probability = " << old_prob << endl;
#endif

      // break if not improving
      if (new_c / old_c0 > 0.995 and tours > NUM_OPTIMIZATION_TOURS)
        break;
    }
    save_coefficients(pr, best_b);
  }
  else
  {
    save_coefficients(pr, b);
  }
}

template <class FT> void Pruner<FT>::optimize_coefficients_prob(/*io*/ vector<double> &pr)
{
  vec b(n), best_b(n);
  FT prob, cost;

  // step 1 call global optimization (without fixing succ. prob)
  optimize_coefficients_evec(pr);
  optimize_coefficients_smooth(pr);
  optimize_coefficients_full(pr);

#ifdef DEBUG_PRUNER_OPTIMIZE
  load_coefficients(b, pr);
  prob = measure_metric(b);
  cost = repeated_enum_cost(b);
  cerr << "# [Stage ho] half-optimization done." << endl;
  cerr << "# [Stage ho]    all_enum_cost = " << cost << endl;
  cerr << "# [Stage fo] single_enum_cost = " << single_enum_cost(b) << endl;
  cerr << "# [Stage ho] succ_probability = " << prob << endl;
  cerr << b << endl;

#endif

  // step 2 achieve target succ. prob
  optimize_coefficients_smooth(pr);
  load_coefficients(b, pr);
  prob = measure_metric(b);
  if (prob <= target)
    optimize_coefficients_prob_incr(pr);
  else
  {
    optimize_coefficients_prob_decr(pr);
  }

  // step 3: some local tweaking
  optimize_coefficients_smooth(pr);
  optimize_coefficients_prob_tune(pr);

#ifdef DEBUG_PRUNER_OPTIMIZE
  load_coefficients(b, pr);
  prob = measure_metric(b);
  cost = repeated_enum_cost(b);
  cerr << "# [Stage fo] full-optimization done." << endl;
  cerr << "# [Stage fo]    all_enum_cost = " << cost << endl;
  cerr << "# [Stage fo] single_enum_cost = " << single_enum_cost(b) << endl;
  cerr << "# [Stage fo] succ_probability = " << prob << endl;
#endif

  save_coefficients(pr, b);
}

/**
 *  optimize either overall cost or cost by fixing probability
 */
template <class FT> void Pruner<FT>::optimize_coefficients(/*io*/ vector<double> &pr)
{
  if (opt_overall)
  {
    optimize_coefficients_cost(pr);
  }
  else
  {
    optimize_coefficients_prob(pr);
  }
}

FPLLL_END_NAMESPACE
