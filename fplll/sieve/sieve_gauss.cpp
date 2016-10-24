
/*
  test small svp challenges
  --------------------------------------
  dim  seed  norm current_status  time
  --------------------------------------
  40    0    1702.46439022965    ok
  50    0    1893.16982862077    ok
  60    0    1943.40088504662    ok
 */

#include "sieve_gauss.h"
#include "sieve_gauss_2sieve.cpp"
#include "sieve_gauss_3sieve.cpp"
#include "sieve_gauss_4sieve.cpp"
#include "wrapper.h"

#define REDUCE_TIMING

/**
 * constructor
 */
template <class ZT, class F>
GaussSieve<ZT, F>::GaussSieve(ZZ_mat<ZT> &B, int alg_arg, bool ver, int seed)
{

  /* stats */
  b             = B;
  nr            = b.get_rows();
  nc            = b.get_cols();
  max_list_size = 0;
  iterations    = 0;
  collisions    = 0;
  reductions    = 0;
  samples       = 0;
  goal_sqr_norm = 0;
  mem_lower     = pow(2.0, 0.18 * nc);
  alg           = alg_arg;
  set_verbose(ver);

  /* sanity check */
  if (alg == 2)
  {
    if (verbose)
      cout << "# [info] running 2-sieve" << endl;
    mult            = 0.1;
    add             = 200.0;
    iterations_step = 200;
  }
  else if (alg == 3)
  {
    if (verbose)
      cout << "# [info] running 3-sieve" << endl;
    mult            = 0.1;
    add             = 100.0;
    iterations_step = 50;
  }
  else if (alg == 4)
  {
    if (verbose)
      cout << "# [info] running 4-sieve" << endl;
    mult            = 0.1;
    add             = 50.0;
    iterations_step = 5;
  }
  else
  {
    cout << " Error, only support 2-, 3- and 4-sieve" << endl;
    exit(1);
  }

  /* clean up list */
  free_list_queue();

  /* initialize sampler */
  Sampler = new KleinSampler<ZT, F>(b, verbose, seed);

  /* initialize list */
  init_list();

  /* further initialization by randomization */
  // init_list_rand();

  /* done initialization */
  max_list_size = List.size();

  /* output stats */
  if (verbose)
  {
    cout << "# [info] done initialization, size(List)=" << List.size() << endl;
    cout << "# [info] done initialization, size(Queue)=" << Queue.size() << endl;
    cout << "# [info] done initialization, mem_est=" << mem_lower << endl;
  }
}

/**
 * deconstructor
 */
template <class ZT, class F> GaussSieve<ZT, F>::~GaussSieve()
{
  free_list_queue();
  free_sampler();
}

/**
 * put matrix vectors to list
 */
template <class ZT, class F> void GaussSieve<ZT, F>::add_mat_list(ZZ_mat<ZT> &B)
{
  Z_NR<ZT> t, current_norm;
  dot_product(best_sqr_norm, B[0], B[0]);
  ListPoint<ZT> *p;

  for (int i = 0; i < nr; ++i)
  {
    p = new_listpoint<ZT>(nc);
    matrix_row_to_list_point(B[i], p);

    // cout << "# [info] init: additing point ";
    // cout << p->v << endl;

    if (alg == 3)
      current_norm = update_p_3reduce(p);
    else if (alg == 2)
      current_norm = update_p_2reduce(p);
    else if (alg == 4)
      current_norm = update_p_4reduce(p);
    else
    {
      cout << " Error, only support 2-, 3- and 4-sieve" << endl;
      exit(1);
    }

    if ((current_norm < best_sqr_norm) && (current_norm > 0))
      // if ((current_norm < best_sqr_norm) )
      best_sqr_norm = current_norm;
  }
}

/**
 * init function (used in constructor)
 */
template <class ZT, class F> void GaussSieve<ZT, F>::init_list() { add_mat_list(b); }

/**
 * init pool of samples (used later)
 */
template <class ZT, class F> void GaussSieve<ZT, F>::init_list_rand()
{
  /* after transformation, the size could be large */
  ZZ_mat<mpz_t> NewZ(nr, nc);
  ZZ_mat<ZT> New(nr, nc);
  mpz_t tmp;
  Z_NR<mpz_t> tmpZ;
  mpz_init(tmp);
  FP_NR<double> c, t;
  Z_NR<ZT> x;
  c = 0.0;
  t = 32.0;

  /* init */
  for (int i = 0; i < nr; i++)
  {
    for (int j = 0; j < nc; j++)
    {
      (b[i][j]).get_mpz(tmp);
      NewZ[i][j] = tmp;
    }
  }

  /* randomization */
  for (int i = 0; i < nr; i++)
  {
    for (int k = 0; k < nr; k++)
    {
      if (i != k)
      {
        x = sample_z_basic_alt<ZT, FP_NR<double>>(c, t);
        x.get_mpz(tmp);
        tmpZ = tmp;
        (NewZ[i]).addmul(NewZ[k], tmpZ, (NewZ[k]).size());
      }
    }
  }

  /* reduce */
  lll_reduction(NewZ, LLL_DEF_DELTA, LLL_DEF_ETA, LM_FAST);

  /* set */
  for (int i = 0; i < nr; i++)
  {
    for (int j = 0; j < nc; j++)
    {
      tmpZ = (NewZ[i][j]).get_data();
      tmpZ.get_mpz(tmp);
      New[i][j] = tmp;
    }
  }

  /* add to list */
  add_mat_list(New);
  mpz_clear(tmp);
}

/**
 * free list and queue
 */
template <class ZT, class F> void GaussSieve<ZT, F>::free_list_queue()
{
  /* clean list */
  typename list<ListPoint<ZT> *>::iterator lp_it;
  for (lp_it = List.begin(); lp_it != List.end(); ++lp_it)
    del_listpoint<ZT>(*lp_it);
  List.clear();

  /* clean queue */
  while (!Queue.empty())
  {
    del_listpoint<ZT>(Queue.front());
    Queue.pop();
  }

  /* clean priority queue */
  while (!Queue_Samples.empty())
  {
    del_listpoint<ZT>(Queue_Samples.top());
    Queue_Samples.pop();
  }
}

/**
 * free sampler
 */
template <class ZT, class F> void GaussSieve<ZT, F>::free_sampler() { delete Sampler; }

/**
 * set targeted norm^2
 */
template <class ZT, class F> void GaussSieve<ZT, F>::set_goal_norm2(Z_NR<ZT> norm)
{
  goal_sqr_norm = norm;
}

/**
 * set verbose
 */
template <class ZT, class F> void GaussSieve<ZT, F>::set_verbose(bool ver) { verbose = ver; }

/**
 * print current info
 */
template <class ZT, class F> void GaussSieve<ZT, F>::print_curr_info()
{
  if (verbose)
  {
    if (iterations % iterations_step == 0)
    {
      cout << "# [info] [" << iterations << "] cols=" << collisions;
      cout << " (" << mult * max_list_size + add << ")";
      cout << " reds=" << reductions;
      cout << " |L|=" << List.size();
      cout << " |Q|=" << Queue.size();
      cout << " |samples|=" << samples;
      cout << " |sv|^2=" << List.front()->norm;
      cout << endl;
      cout << std::flush;
    }
  }
}

/**
 * print final info
 */
template <class ZT, class F> void GaussSieve<ZT, F>::print_final_info()
{
  long first_size            = 0;
  vector<long>::iterator it2 = iters_ls.begin();
  for (typename vector<Z_NR<ZT>>::iterator it1 = iters_norm.begin(); it1 != iters_norm.end();
       ++it1, ++it2)
  {
    if ((*it1) == best_sqr_norm)
    {
      first_size = (*it2);
      break;
    }
  }
  if (verbose)
  {
    cout << "# [****] done!" << endl;
    cout << "# [info] [" << iterations << "] cols=" << collisions;
    cout << " (" << mult * max_list_size + add << ")";
    cout << " reds=" << reductions;
    cout << " |L|=" << List.size();
    cout << " |Q|=" << Queue.size();
    cout << " |samples|=" << samples << endl;
    cout << "# [info] max(|L|)=" << max_list_size;
    cout << " log2(max|L|)/n=" << log2(max_list_size) / nc << endl;
    cout << "# [info] true max|L| = " << first_size << endl;
    ;
    cout << "# [info] true log2(max|L|)/n = " << log2(first_size) / nc << endl;
    cout << "# [info] sv is" << endl;
  }
  cout << List.front()->v << endl;
  if (verbose)
  {
    final_norm.set_z(best_sqr_norm);
    final_norm.sqrt(final_norm, GMP_RNDN);
    cout << "# [info] |sv| = " << final_norm << " (" << best_sqr_norm << ")" << endl;
  }
}

template <class ZT, class F> NumVect<Z_NR<ZT>> GaussSieve<ZT, F>::return_first()
{
  return List.front()->v;
}

template class GaussSieve<mpz_t, FP_NR<double>>;
template class GaussSieve<long, FP_NR<double>>;
