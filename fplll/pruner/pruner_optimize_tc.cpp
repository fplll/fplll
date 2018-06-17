#include "fplll.h"

FPLLL_BEGIN_NAMESPACE

#define BALANCE_HEURISTIC_PRUNER_OPTIMIZE
//#define DEBUG_PRUNER_OPTIMIZE_TC

/**
 *  optimize with constrains that b_i = b_{i+1} for even i.
 */
template <class FT> void Pruner<FT>::optimize_coefficients_evec(/*io*/ vector<double> &pr)
{
  evec b(d);

  // load coefficients
  if (flags & PRUNER_START_FROM_INPUT)
  {
    load_coefficients(b, pr);
  }

  // greedy method
  if (!(flags & PRUNER_START_FROM_INPUT))
  {
    greedy(b);
#ifdef DEBUG_PRUNER_OPTIMIZE_TC
    cerr << "# [Greedy]" << endl;
    cerr << b << endl;
    cerr << "# [Greedy] single_enum_cost  = " << single_enum_cost(b) << endl;
    cerr << "# [Greedy] succ_probability  = " << measure_metric(b) << endl;
    cerr << "# [Greedy]    all_enum_cost  = " << repeated_enum_cost(b) << endl;
#endif
  }

  // greedy method for min pruning coefficients.
  if (flags & (PRUNER_GRADIENT | PRUNER_NELDER_MEAD))
  {
    preproc_cost *= .1;
    greedy(min_pruning_coefficients);

    // The aim is to get a lower bound for the pruning parameter
    // note this lower bound will be used in the enforce()
    // In the case of fixed input probability, check whether this
    // min_pruning_coefficiens is small enough. Otherwise, reduce
    // it further. This is important since otherwise one may never
    // achieve the target probability.
    if (!opt_overall)
    {
      vector<double> pr(n);
      save_coefficients(pr, min_pruning_coefficients);
      if (measure_metric(min_pruning_coefficients) > target)
      {
        fill(min_pruning_coefficients.begin(), min_pruning_coefficients.end(), 0.);
        optimize_coefficients_decr_prob(pr);
      }
      load_coefficients(min_pruning_coefficients, pr);
    }
    preproc_cost *= 10;
  }

  // 2. gradient method // modify this to becomes an independent method!!!!
  if (flags & PRUNER_GRADIENT)
  {
    if (verbosity)
    {
      cerr << "\nGradient descent start (dim=" << n << ")" << endl;
    }
    gradient_descent(b);
#ifdef DEBUG_PRUNER_OPTIMIZE_TC
    cerr << "# [Descent]" << endl;
    cerr << b << endl;
    cerr << "# [Descent] single_enum_cost  = " << single_enum_cost(b) << endl;
    cerr << "# [Descent] succ_probability  = " << measure_metric(b) << endl;
    cerr << "# [Descent]    all_enum_cost  = " << repeated_enum_cost(b) << endl;
#endif
  };

  if (flags & PRUNER_NELDER_MEAD)
  {
    if (verbosity)
    {
      cerr << "\nNelder-Mead start (dim=" << n << ")" << endl;
    }
    while (nelder_mead_step(b))
    {
    };
#ifdef DEBUG_PRUNER_OPTIMIZE_TC
    cerr << "# [Nelder-Mead]" << endl;
    cerr << b << endl;
    cerr << "# [Nelder-Mead] single_enum_cost  = " << single_enum_cost(b) << endl;
    cerr << "# [Nelder-Mead] succ_probability  = " << measure_metric(b) << endl;
    cerr << "# [Nelder-Mead]    all_enum_cost  = " << repeated_enum_cost(b) << endl;
#endif
  };
  save_coefficients(pr, b);
}

/**
 *  optimize without constrains b_i = b_{i+1} for even i.
 */
template <class FT> void Pruner<FT>::optimize_coefficients_full(/*io*/ vector<double> &pr)
{
  vec b(n);

  // always load coefficients since this is used after optimize_coefficients_evec
  load_coefficients(b, pr);

  // gradient method
  if (flags & PRUNER_GRADIENT)
  {
    if (verbosity)
    {
      cerr << "\nGradient descent start (dim=" << n << ")" << endl;
    }

    gradient_descent(b);
#ifdef DEBUG_PRUNER_OPTIMIZE_TC
    cerr << "# [Descent]" << endl;
    cerr << b << endl;
    cerr << "# [Descent] single_enum_cost  = " << single_enum_cost(b) << endl;
    cerr << "# [Descent] succ_probability  = " << measure_metric(b) << endl;
    cerr << "# [Descent]    all_enum_cost  = " << repeated_enum_cost(b) << endl;
#endif
  };

  // Nelder-Mead method
  if (flags & PRUNER_NELDER_MEAD)
  {
    if (verbosity)
    {
      cerr << "\nNelder-Mead start (dim=" << n << ")" << endl;
    }
    while (nelder_mead_step(b))
    {
    };

#ifdef DEBUG_PRUNER_OPTIMIZE_TC
    cerr << "# [Nelder-Mead]" << endl;
    cerr << b << endl;
    cerr << "# [Nelder-Mead] single_enum_cost  = " << single_enum_cost(b) << endl;
    cerr << "# [Nelder-Mead] succ_probability  = " << measure_metric(b) << endl;
    cerr << "# [Nelder-Mead]    all_enum_cost  = " << repeated_enum_cost(b) << endl;
#endif
  };

  save_coefficients(pr, b);
}

/**
 * Tweaking of pruning coefficients in neighborhood: try to reduce
 * single enum cost (hopefull reduce the overall running time)
 */
template <class FT>
void Pruner<FT>::optimize_coefficients_local_adjust_decr_single(/*io*/ vector<double> &pr)
{
  int maxi, lasti, consecutive_fails;
  double improved_ratio, current_max = 0.0;
  FT old_cf, old_cfs, new_cf, old_b;
  vector<double> detailed_cost(n);
  vector<double> slices(n, 10.0);  // (b[i+1] - b[i])/slice will be used as step
  vector<int> thresholds(n, 3);
  vec b(n);

  load_coefficients(b, pr);

  lasti = -1;             // last failed index, make sure we do not try it again in
                          // the next time
  consecutive_fails = 0;  // number of consecutive failes; break if
                          // reaches it

  improved_ratio = 0.995;  // if reduced by 0.995, descent

  while (1)
  {

    // old cost
    old_cf = target_function(b);

    // find bottleneck index
    old_cfs = single_enum_cost(b, &(detailed_cost));

// heuristic
#ifdef BALANCE_HEURISTIC_PRUNER_OPTIMIZE
    if (old_cfs < sqrt(old_cf) / 10.0)
      break;
#endif

    current_max = 0.0;
    maxi        = 0;
    for (int i = 0; i < n; i++)
    {
      if ((i != (n - lasti - 1)) && (thresholds[n - i - 1] > 0))
      {
        if (detailed_cost[i] > current_max)
        {
          current_max = detailed_cost[i];
          maxi        = i;
        }
      }
    }

    // b[ind] is the one to be reduced
    int ind = n - maxi - 1;
    old_b   = b[ind];
    if (ind != 0)
    {
      b[ind] = b[ind] - (b[ind] - b[ind - 1]) / slices[ind];
    }
    else
    {
      break;
    }

    // new cost
    new_cf = target_function(b);

    // if not improved -- recover
    if (new_cf >= (old_cf * improved_ratio))
    {
      b[ind] = old_b;
      lasti  = ind;
      thresholds[lasti]--;
      consecutive_fails++;
    }
    else
    {
      // cerr << " improved from " << old_cf << " to  " << new_cf << endl;
      if (slices[ind] < 1024)
        slices[ind]     = slices[ind] * 1.05;
      consecutive_fails = 0;
    }

    // quit after 10 consecutive failes
    if (consecutive_fails > 10)
    {
      break;
    }
  }

#ifdef DEBUG_PRUNER_OPTIMIZE_TC
  cerr << "# [TuningCost]" << endl;
  cerr << b << endl;
  cerr << "# [TuningCost] all_enum_cost    = " << repeated_enum_cost(b) << endl;
  cerr << "# [TuningCost] succ_probability = " << measure_metric(b) << endl;
#endif

  save_coefficients(pr, b);
}

/**
 * Tuning of pruning coefficients: try to increase prob by increasing
 * those coefficients b[i] inversely-weighted by their level cost.
 */
template <class FT>
void Pruner<FT>::optimize_coefficients_local_adjust_incr_prob(/*io*/ vector<double> &pr)
{
  int trials, tours, maxi, ind;
  FT old_cf, old_cf0, old_cfs, new_cf, old_b;
  double current_max;
  vector<double> detailed_cost(n);
  vector<double> slices(n, 10.0);  // (b[i+1] - b[i])/slice will be used as step
  vec b(n);
  load_coefficients(b, pr);

  // initial cost
  old_cf0 = target_function(b);

  tours = 0;
  while (1)
  {

    tours++;

    // old cost
    old_cf = target_function(b);
    // find bottleneck index
    old_cfs     = single_enum_cost(b, &(detailed_cost));
    current_max = 0.0;
    maxi        = 0;
    for (int i = 0; i < n; i++)
    {
      if (detailed_cost[i] > current_max)
      {
        current_max = detailed_cost[i];
        maxi        = i;
      }
    }
    ind = n - maxi - 1;
    if (ind <= 1)
      break;

#ifdef BALANCE_HEURISTIC_PRUNER_OPTIMIZE
    if (old_cfs > sqrt(old_cf) / 10.0)
      break;
#endif

    for (int i = ind; i >= 1; --i)
    {

      if (b[i] <= b[i - 1])
        continue;

      trials = 0;

      // fixed i-1, trying to increase b[i-1]
      while (1)
      {
        // old cost
        old_cf = target_function(b);

        // try increase
        old_b    = b[i - 1];
        b[i - 1] = b[i - 1] + (b[i] - b[i - 1]) / slices[i - 1];

        // new cost
        new_cf = target_function(b);

        // cerr << " i = " << i << " old_cf = " << old_cf << " new_cf = " <<
        // new_cf << endl;

        // if not improved -- recover
        if (new_cf >= (old_cf * 1.2))
        {
          b[i - 1] = old_b;
          break;
        }
        else
        {
          if (slices[i - 1] < 1024)
            slices[i - 1] = slices[i - 1] * 1.2;
        }
        trials++;
        if (trials >= 10)
          break;
      }
    }

    new_cf = target_function(b);
    if (new_cf > (old_cf0 * 1.1) || tours > 4)
      break;
  }

#ifdef DEBUG_PRUNER_OPTIMIZE_TC
  cerr << "# [TuningProb]" << endl;
  cerr << b << endl;
  cerr << "# [TuningProb] all_enum_cost    = " << repeated_enum_cost(b) << endl;
  cerr << "# [TuningProb] succ_probability = " << measure_metric(b) << endl;
#endif

  save_coefficients(pr, b);
}

/**
 * try to smooth the curve if there are obvious discountinuties between
 * consecutive indices
 */
template <class FT>
void Pruner<FT>::optimize_coefficients_local_adjust_smooth(/*io*/ vector<double> &pr)
{
  vec b(n);
  FT lr, rr;
  FT th = 1.0 / n;
  load_coefficients(b, pr);

  for (int i = 1; i < n - 1; ++i)
  {

    lr = b[i] / b[i - 1];
    rr = b[i + 1] / b[i];

    if ((rr / lr > 1.25) || (rr / lr < 0.8))
    {
      b[i] = sqrt(b[i - 1] * b[i + 1]);
    }

    if ((b[i + 1] - b[i]) > th || (b[i] - b[i - 1]) > th)
    {
      b[i] = (b[i - 1] + b[i + 1]) / 2.0;
    }
  }

  save_coefficients(pr, b);
}

/**
 * Optimize b with the greedy method
 */
template <class FT> void Pruner<FT>::greedy(evec &b)
{
  // Do not call enforce in this function, as min_pruning_bounds may not have been set
  // Indeed, the min_pruning_bound should now based on greedy.
  if (!shape_loaded)
  {
    throw std::invalid_argument("Error: No basis shape was loaded");
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

/**
 * Optimize b with the gradient descent
 */
template <class FT> int Pruner<FT>::gradient_descent(/*io*/ vec &b)
{

  FT old_epsilon  = epsilon;
  FT old_min_step = min_step;
  int trials      = 0;

  while (1)
  {
    int ret = gradient_descent_step(b);
    if (ret == 0)
      break;
    else if (ret < 0)
    {
      epsilon  = epsilon * 0.9;
      min_step = min_step * 0.9;
      trials++;
      if (trials >= 5)
        break;
    }
    else
    {
      trials--;
      continue;
    }
  }

  epsilon  = old_epsilon;
  min_step = old_min_step;
  return 0;
}

/**
 * One gradient descent step
 */
template <class FT> int Pruner<FT>::gradient_descent_step(/*io*/ vec &b)
{
  int dn    = b.size();
  FT cf     = target_function(b);
  FT old_cf = cf;
  vec new_b(dn);
  vector<double> pr(dn);
  vec gradient(dn);
  target_function_gradient(b, gradient);
  FT norm = 0.0;

  // normalize the gradient
  for (int i = 0; i < dn; ++i)
  {
    norm += gradient[i] * gradient[i];
    new_b[i] = b[i];
  }

  if (verbosity)
  {
    cerr << "  Gradient descent step starts at cf=" << cf << endl;
  }

  norm /= (double)dn;
  norm = sqrt(norm);

  if (verbosity)
  {
    cerr << "  Gradient norm " << norm << endl;
  }

  if (norm <= 0.)
    return 0;

  for (int i = 0; i < dn; ++i)
  {
    gradient[i] /= norm;
  }
  FT new_cf;

  FT step = min_step;
  int j;

  for (j = 0;; ++j)
  {
    if (step > dn)
    {
      return -1;
      // throw std::runtime_error("Infinite loop in pruner gradient_descent_step");
    }

    for (int i = 0; i < dn; ++i)
    {
      new_b[i] = new_b[i] + step * gradient[i];
    }

    enforce(new_b);

    new_cf = target_function(new_b);

    if (new_cf >= cf)
    {
      break;
    }
    b  = new_b;
    cf = new_cf;
    step *= step_factor;
  }

  if (verbosity)
  {
    cerr << "  Gradient descent step ends after " << j << " mini-steps at cf=" << cf << endl;
  }

  if (cf > old_cf * min_cf_decrease)
  {
    return 0;
  }
  return j;
}

/**
 * Optimize b with Nelder-Mead method. Following the notation of
 *    https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
 */
#define ND_ALPHA 1
#define ND_GAMMA 2
#define ND_RHO 0.5
#define ND_SIGMA 0.5
#define ND_INIT_WIDTH 0.01

template <class FT> int Pruner<FT>::nelder_mead_step(/*io*/ vec &b)
{
  int dn = b.size();
  int l  = dn + 1;
  evec tmp_constructor(dn);
  vector<evec> bs(l);  // The simplexe (d+1) vector of dim d
  FT *fs = new FT[l];  // Values of f at the simplex vertices

  for (int i = 0; i < l; ++i)  // Intialize the simplex
  {
    bs[i] = b;  // Start from b
    if (i < dn)
    {
      bs[i][i] += (bs[i][i] < .5) ? ND_INIT_WIDTH : -ND_INIT_WIDTH;
    }
    enforce(bs[i]);
    fs[i] = target_function(bs[i]);  // initialize the value
  }

  FT init_cf = fs[l - 1];

  vec bo(dn);  // centeroid
  FT fo;       // value at the centroid

  FT fs_maxi_last = fs[0];  // value of the last centroid

  if (verbosity)
  {
    cerr << "  Starting nelder_mead cf = " << init_cf << " proba = " << measure_metric(b) << endl;
  }
  unsigned int counter = 0;
  int mini = 0, maxi = 0, maxi2 = 0;
  while (1)  // Main loop
  {
    mini = maxi = maxi2 = 0;
    for (int i = 0; i < dn; ++i)
      bo[i]    = bs[0][i];
    ////////////////
    // step 1. and 2. : Order and centroid
    ////////////////
    for (int i = 1; i < l; ++i)  // determine min and max, and centroid
    {
      mini = (fs[i] < fs[mini]) ? i : mini;
      maxi = (fs[i] > fs[mini]) ? i : maxi;
      for (int j = 0; j < dn; ++j)
      {
        bo[j] += bs[i][j];
      }
    }
    FT tmp;
    tmp = l;
    for (int i = 0; i < dn; ++i)
      bo[i] /= tmp;  // Centroid calculated

    fs_maxi_last = (!counter) ? fs[maxi] : fs_maxi_last;

    // determine min and max, and centroid
    maxi2 += (!maxi);
    for (int i = 1; i < l; ++i)
    {
      maxi2 = ((fs[i] > fs[maxi2]) && (i != maxi)) ? i : maxi2;
    }

    if (enforce(bo))  // Maybe want to skip this test, that may be kinda costly.
    {
      throw std::runtime_error("Concavity says that should not happen.");
    }

    if (verbosity)  // Not sure such verbosity is now of any use
    {
      cerr << "  melder_mead step " << counter << "cf = " << fs[mini]
           << " proba = " << measure_metric(bs[mini]) << " cost = " << single_enum_cost(bs[mini])
           << endl;
      for (int i = 0; i < dn; ++i)
      {
        cerr << ceil(bs[mini][i].get_d() * 1000) << " ";
      }
      cerr << endl;
    }

    ////////////////
    // Stopping condition (Not documented on wikipedia, improvising)
    // I'm not satistifed by it anymore. Exploration needed.
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

    vec br(dn);  // reflected point
    FT fr;       // Value at the reflexion point
    for (int i = 0; i < dn; ++i)
      br[i]    = bo[i] + ND_ALPHA * (bo[i] - bs[maxi][i]);
    enforce(br);
    fr = target_function(br);
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
      vec be(dn);
      FT fe;
      for (int i = 0; i < dn; ++i)
        be[i]    = bo[i] + ND_GAMMA * (br[i] - bo[i]);
      enforce(be);
      fe = target_function(be);
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

    vec bc(dn);
    FT fc;
    for (int i = 0; i < dn; ++i)
      bc[i]    = bo[i] + ND_RHO * (bs[maxi][i] - bo[i]);
    enforce(bc);
    fc = target_function(bc);
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
      for (int i = 0; i < dn; ++i)
      {
        bs[j][i] = bs[mini][i] + ND_SIGMA * (bs[j][i] - bs[mini][i]);
      }
      enforce(bs[j]);
      fs[j] = target_function(bs[j]);  // initialize the value
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

FPLLL_END_NAMESPACE
