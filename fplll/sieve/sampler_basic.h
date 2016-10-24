/* "Trapdoors for Hard Lattices and New Cryptographic Constructions"
   Craig Gentry, Chris Peikert, Vinod Vaikuntanathan */

#ifndef FPLLL_SIEVE_SAMPLER_BASIC_H
#define FPLLL_SIEVE_SAMPLER_BASIC_H

#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"

using namespace std;
using namespace fplll;

template <class ZT, class F> class KleinSampler
{

public:
  bool verbose;

  KleinSampler(ZZ_mat<ZT> &B, bool verbose, int seed);
  ~KleinSampler();
  void print_param();
  Z_NR<ZT> sample_z(F c, F s);
  Z_NR<ZT> sample_z_basic(F c, F s);
  void set_verbose(bool verbose);

  NumVect<Z_NR<ZT>> sample();

private:
  /**
   * GSO object of B
   */
  int nr, nc;
  double logn2;
  ZZ_mat<ZT> b, u, u_inv;
  MatGSO<Z_NR<ZT>, F> *pGSO;
  Matrix<Z_NR<ZT>> g;
  Matrix<F> mu, r;

  /* input parameter s */
  F maxbistar2;

  /* input parameter s */
  F s2;

  /* input parameter s */
  F c;

  /* ratio t = logn */
  F t;

  /* variances */
  NumVect<F> *s_prime;
};

/**
 * sampling Z by rejection sampling (not a class member)
 */
template <class ZT, class F> Z_NR<ZT> sample_z_basic_alt(F c, F s)
{
  F min, max, st, range, tmp, tmp1;
  double r, e;

  /* (c \pm s*t) for t \approx logn */
  st = s;
  F t;
  t = 2.0;
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

#endif
