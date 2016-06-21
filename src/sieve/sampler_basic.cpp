#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sampler_basic.h"


/**************************
 *  class KleinSampler
 **************************/


template<class ZT, class F>
KleinSampler<ZT, F>::KleinSampler (ZZ_mat<ZT> &B)
{
  /* set dimensions */
  b = B;
  nr = b.getRows();
  nc = b.getCols();
  //t = log(nr);
  t = 2;
  logn2 = log(nr)*log(nr);

  /* gso, flag 1 to have g matrix valid */
  pGSO = new MatGSO<Z_NR<ZT>, F> (b, u, uInv, 1);

  pGSO->updateGSO();
  mu = pGSO->getMuMatrix();
  r = pGSO->getRMatrix();
  g = pGSO->getGMatrix();
  
  /* compute variances for sampling */
  maxbistar2 = pGSO->getMaxBstar();
  s2.mul_d (maxbistar2, logn2, GMP_RNDN);
  s_prime = new NumVect<F> (nr);
  F tmp;
  for (int i = 0; i < nr; i++) {
    tmp.set_z(g(i,i));
    ((*s_prime)[i]).div(maxbistar2, tmp);
    ((*s_prime)[i]).sqrt((*s_prime)[i], GMP_RNDN);
  }

  /* verbose */
  print_param();
  gmp_randinit_default(state);
}


template<class ZT, class F>
KleinSampler<ZT, F>::~KleinSampler ()
{
  gmp_randclear(state);
  delete pGSO;
  delete s_prime;
}


template<class ZT, class F>
void KleinSampler<ZT, F>::print_param ()
{
  /*
    cout << b << endl;
    cout << r << endl;
    cout << g << endl;
  */
  cout << "# [info] nc = " << nc << endl;
  cout << "# [info] nr = " << nr << endl;
  cout << "# [info] t = log(nr) = " << t << endl;
  cout << "# [info] maxbistar2 = " << maxbistar2 << endl;
  /*
    for (int i = 0; i < nr; i++) {
    cout << "# [info]  B_" << i << "^* = " << g(i,i) << 
    ", (s')^2 = " << ((*s_prime)[i]) << endl;
    }
  */
}


/**
 * sampling Z by rejection sampling
 */
template<class ZT, class F>
Z_NR<ZT> KleinSampler<ZT, F>::sampleZ_basic (F c, F s)
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

  /*
  cout << "st is " << st << endl;
  cout << "min is " << min << endl;
  cout << "max is " << max << endl;
  */

  Z_NR<ZT> x;
  while(1) {
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
 * (long, double) use function from dgs
 */
template<>
Z_NR<long> KleinSampler<long, FP_NR<double> >::sampleZ_dgs
( FP_NR<double> c,
  FP_NR<double> s,
  long flag )
{
  Z_NR<long> rz;
  double s_ = s.get_d();
  double c_ = c.get_d();
  dgs_disc_gauss_alg_t alg = DGS_DISC_GAUSS_UNIFORM_TABLE;
  dgs_disc_gauss_dp_t *gen = dgs_disc_gauss_dp_init(s_, c_, 6, alg);
  rz = gen->call(gen);
  dgs_disc_gauss_dp_clear(gen);
  return rz;
}


/**
 * (mpz_t, mpfr_t) use function from dgs
 */
template<>
Z_NR<mpz_t> KleinSampler<mpz_t, FP_NR<mpfr_t> >::sampleZ_dgs
( FP_NR<mpfr_t> c,
  FP_NR<mpfr_t> s,
  mpz_t flag )
{
  mpfr_set_default_prec(80);
  Z_NR<mpz_t> rz;
  mpfr_t s_, c_;
  mpz_t r;
  mpz_init(r);
  mpfr_init(s_);
  mpfr_init(c_);
  s.get_mpfr(s_, GMP_RNDN);
  c.get_mpfr(c_, GMP_RNDN);
  dgs_disc_gauss_alg_t alg = DGS_DISC_GAUSS_UNIFORM_TABLE;
  dgs_disc_gauss_mp_t *gen = dgs_disc_gauss_mp_init(s_, c_, 6, alg);
  gen->call(r, gen, state);
  rz = r;
  dgs_disc_gauss_mp_clear(gen);
  mpz_clear(r);
  mpfr_clear(s_);
  mpfr_clear(c_);
  return rz;
}


/**
 * (mpz_t, double) use sampler_basic
 */
template<>
Z_NR<mpz_t> KleinSampler<mpz_t, FP_NR<double> >::sampleZ_dgs
( FP_NR<double> c,
  FP_NR<double> s,
  mpz_t flag )
{
  return sampleZ_basic (c, s);
}


/**
 * support three modes:
 *   long, double
 *   mpz_t, double
 *   mpz_t, mpfr_t
 */
template<class ZT, class F>
Z_NR<ZT> KleinSampler<ZT, F>::sampleZ (F c, F s)
{
  ZT flag;
  flag = 0;
  return sampleZ_dgs (c, s, flag);
}


template<class ZT, class F>
NumVect<Z_NR<ZT> > KleinSampler<ZT, F>::sample ()
{
  NumVect<Z_NR<ZT> > vec(nr);
  NumVect<F> ci(nr);
  F tmp;
  Z_NR<ZT> tmpz;

  // while(1) {
    
    for (int i = 0; i < nr; i++) {
      ci[i] = 0;
      vec[i] = 0;
    }

    for (int i = nr - 1; i >= 0; i--) {
      tmpz = sampleZ (ci[i], (*s_prime)[i]);
      (ci[i]).set_z(tmpz);
      for (int j = 0; j < i; j++) {
        tmp.mul(ci[i], mu(i,j), GMP_RNDN);
        (ci[j]).sub(ci[j], tmp, GMP_RNDN);
      }
    }

    //lp->norm = 0;
    for (int i = 0; i < nc; i++) {
      for (int j = 0; j < nr; j++) {
        tmpz.set_f(ci[j]);
        tmpz.mul(tmpz, b(j,i));
        (vec[i]).add(vec[i], tmpz);
      }
      //lp->norm += lp->v[i] * lp->v[i];
    }
/*
    if (!vec.is_zero())
      break;
  }
*/
  return vec;
}


template class KleinSampler<long, FP_NR<double> >;
template class KleinSampler<mpz_t, FP_NR<double> >;
template class KleinSampler<mpz_t, FP_NR<mpfr_t> >;

