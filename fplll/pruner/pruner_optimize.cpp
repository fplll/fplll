#include "fplll.h"

FPLLL_BEGIN_NAMESPACE

#define DEBUG_PRUNER_OPTIMIZE
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
  cout << "# [Stage ho] half-optimization done." << endl;
  cout << "# [Stage ho]  all_enum_cost    is " << old_c0 << endl;
  cout << "# [Stage ho]  single_enum_cost is " << old_sc << endl;
  cout << "# [Stage ho]         succ_prob is " << old_prob << endl;
#endif

  // step 2 use full coefficients
  int tours = 0;
  while (1)
  {

    tours++;

    // initial cost
    load_coefficients(b, pr);
    old_c0 = repeated_enum_cost(b);

    // step 2.1 full tuning
    optimize_coefficients_tune_cost(pr);
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
    cout << "# [Stage ft] full-tuning done." << endl;
    cout << "# [Stage ft]  all_enum_cost    is " << old_c1 << endl;
    cout << "# [Stage ft]  single_enum_cost is " << old_sc << endl;
    cout << "# [Stage ft]         succ_prob is " << old_prob << endl;
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
    cout << "# [Stage fio] final optimization done." << endl;
    cout << "# [Stage fio]  all_enum_cost    is " << new_c << endl;
    cout << "# [Stage fio]  single_enum_cost is " << old_sc << endl;
    cout << "# [Stage fio]         succ_prob is " << old_prob << endl;
#endif

    if (new_c / old_c0 > 0.995 and tours > NUM_OPTIMIZATION_TOURS)
      break;
  }

  cout << best_b << endl;

  // done
  save_coefficients(pr, best_b);
}

template <class FT> void Pruner<FT>::optimize_coefficients_prob(/*io*/ vector<double> &pr)
{
  vec b(n), best_b(n);
  FT prob, cost;

  // step 1 use half coefficients only
  optimize_coefficients_evec(pr);
  optimize_coefficients_smooth(pr);
  optimize_coefficients_full(pr);

#if 1
  load_coefficients(b, pr);
  prob = measure_metric(b);
  cost = repeated_enum_cost(b);
  cerr << "# [Stage ho] half-optimization done." << endl;
  cerr << "# [Stage ho]  all_enum_cost    is " << cost << endl;
  cerr << "# [Stage ho]  svp_probability is " << prob << endl;
#endif
  // step 2 use full coefficients
  optimize_coefficients_smooth(pr);

  load_coefficients(b, pr);
  prob = measure_metric(b);

  if (prob <= target)
    optimize_coefficients_prob_incr(pr);
  else
  {
    optimize_coefficients_prob_decr(pr);
  }

  // step 3: adjust
  optimize_coefficients_smooth(pr);
  optimize_coefficients_prob_tune(pr);

#if 1
  load_coefficients(b, pr);
  prob = measure_metric(b);
  cost = repeated_enum_cost(b);
  cerr << "# [stage] half-optimization done." << endl;
  cerr << "# [stage]  all_enum_cost    is " << cost << endl;
  cerr << "# [stage]  svp_probability is " << prob << endl;
#endif
}

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
