#include "sampler_basic.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/**************************
 *  class KleinSampler
 **************************/

template <class ZT, class F> KleinSampler<ZT, F>::KleinSampler(ZZ_mat<ZT> &B, bool ver, int seed)
{
  /* set dimensions */
  b  = B;
  nr = b.get_rows();
  nc = b.get_cols();
  // t = log(nr);
  t     = 2;
  logn2 = log(nr) * log(nr);

  /* gso, flag 1 to have g matrix valid */
  pGSO = new MatGSO<Z_NR<ZT>, F>(b, u, u_inv, 1);

  pGSO->update_gso();
  mu = pGSO->get_mu_matrix();
  r  = pGSO->get_r_matrix();
  g  = pGSO->get_g_matrix();

  /* compute variances for sampling */
  maxbistar2 = pGSO->get_max_bstar();

  s2.mul_d(maxbistar2, logn2, GMP_RNDN);
  s_prime = new NumVect<F>(nr);

  F tmp;
  for (int i = 0; i < nr; i++)
  {
    tmp.set_z(g(i, i));
    ((*s_prime)[i]).div(maxbistar2, tmp);
    ((*s_prime)[i]).sqrt((*s_prime)[i], GMP_RNDN);
  }

  /* verbose */
  srand(seed);
  set_verbose(ver);
  print_param();
}

template <class ZT, class F> KleinSampler<ZT, F>::~KleinSampler()
{
  delete pGSO;
  delete s_prime;
}

/**
 * set verbose
 */
template <class ZT, class F> void KleinSampler<ZT, F>::set_verbose(bool ver) { verbose = ver; }

template <class ZT, class F> void KleinSampler<ZT, F>::print_param()
{
  if (verbose)
  {
    cout << "# [info] nc = " << nc << endl;
    cout << "# [info] nr = " << nr << endl;
    cout << "# [info] t = log(nr) = " << t << endl;
    cout << "# [info] maxbistar2 = " << maxbistar2 << endl;
  }
}

/**
 * sampling Z by rejection sampling
 */
template <class ZT, class F> Z_NR<ZT> KleinSampler<ZT, F>::sample_z_basic(F c, F s)
{
  F min, max, st, range, tmp, tmp1;
  double r, e;

  /* (c \pm s*t) for t \approx logn */
  st = s;
  st.mul(st, t, GMP_RNDN);
  min.sub(c, st, GMP_RNDN);
  max.add(c, st, GMP_RNDN);
  min.rnd(min);
  max.rnd(max);
  range.sub(max, min, GMP_RNDN);

  Z_NR<ZT> x;
  while (1)
  {
    r = double(rand()) / RAND_MAX;
    tmp.mul_d(range, r, GMP_RNDN);
    tmp.rnd(tmp);
    tmp.add(tmp, min, GMP_RNDN);
    x.set_f(tmp);
    tmp1.sub(tmp, c, GMP_RNDN);
    tmp1.mul(tmp1, tmp1, GMP_RNDN);
    tmp1.mul_d(tmp1, -M_PI, GMP_RNDN);
    tmp.mul(s, s, GMP_RNDN);
    tmp1.div(tmp1, tmp, GMP_RNDN);
    e = tmp1.get_d(GMP_RNDN);
    r = exp(e);
    if ((double(rand()) / RAND_MAX) <= r)
      return x;
  }
}

/**
 * support three modes:
 *   long, double
 *   mpz_t, double
 *   mpz_t, mpfr_t
 */
template <class ZT, class F> Z_NR<ZT> KleinSampler<ZT, F>::sample_z(F c, F s)
{
  return sample_z_basic(c, s);
}

template <class ZT, class F> NumVect<Z_NR<ZT>> KleinSampler<ZT, F>::sample()
{
  NumVect<Z_NR<ZT>> vec(nc);
  NumVect<F> ci(nc);
  F tmp;
  Z_NR<ZT> tmpz;

  for (int i = 0; i < nc; i++)
    vec[i]   = 0.0;
  for (int i = 0; i < nr; i++)
    ci[i]    = 0.0;

  for (int i = nr - 1; i >= 0; i--)
  {
    tmpz = sample_z(ci[i], (*s_prime)[i]);
    (ci[i]).set_z(tmpz);
    for (int j = 0; j < i; j++)
    {
      tmp.mul(ci[i], mu(i, j), GMP_RNDN);
      (ci[j]).sub(ci[j], tmp, GMP_RNDN);
    }
  }

  // lp->norm = 0;
  for (int i = 0; i < nc; i++)
  {
    for (int j = 0; j < nr; j++)
    {
      tmpz.set_f(ci[j]);
      tmpz.mul(tmpz, b(j, i));
      (vec[i]).add(vec[i], tmpz);
    }
    // lp->norm += lp->v[i] * lp->v[i];
  }
  return vec;
}

template class KleinSampler<long, FP_NR<double>>;
template class KleinSampler<mpz_t, FP_NR<double>>;
template class KleinSampler<mpz_t, FP_NR<mpfr_t>>;
