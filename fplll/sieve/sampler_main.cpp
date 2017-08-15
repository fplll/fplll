#include "sampler_basic.h"

int main()
{
  srand(time(NULL));
  int dim = 40;

#if 1
  /* long matrix, double */
  ZZ_mat<mpz_t> M1_Z(dim, dim);
  M1_Z.gen_GM(20 * dim);
  lllReduction(M1_Z, 0.99, 0.51, LM_WRAPPER);
  ZZ_mat<long> M1(dim, dim);
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      M1(i, j) = M1_Z(i, j).get_si();
  FP_NR<double> c, t;
  Z_NR<long> x1;
  long flag = 0;
  c         = 0.0;
  t         = 1000.0;
  KleinSampler<long, FP_NR<double>> S1(M1);

  /* sample Z */
  for (int i = 0; i < 1000; i++)
  {
    x1 = S1.sample_z_dgs(c, t, flag);
    // cout << x1 << endl;
  }
  cout << "# [done] (long, double) sample_z_dgs()" << endl;

  /* sample vec */
  NumVect<Z_NR<long>> vec;
  clock_t stime = clock();
  for (int i = 0; i < 10000; i++)
  {
    vec = S1.sample();
    // cout << (vec) << endl;
  }
  cout << "# [done] (long, double) sample()" << endl;
  clock_t etime = clock();
  double secs   = (etime - stime) / (double)CLOCKS_PER_SEC;
  cout << "# [info] (long double) took time " << secs << " s" << endl;
#endif

#if 1
  /* mpz_t matrix, double */
  ZZ_mat<mpz_t> M2(dim, dim);
  M2.gen_GM(20 * dim);
  lllReduction(M2, 0.99, 0.51, LM_WRAPPER);
  Z_NR<mpz_t> x2;
  KleinSampler<mpz_t, FP_NR<double>> S2(M2);

  /* sample Z */
  for (int i = 0; i < 1000; i++)
  {
    x2 = S2.sample_z(c, t);
    // cout << x2 << endl;
  }
  cout << "# [done] (mpz_t, double) sample_z_basic()" << endl;

  /* sample vec */
  NumVect<Z_NR<mpz_t>> vec2;
  stime = clock();
  for (int i = 0; i < 10000; i++)
  {
    vec2 = S2.sample();
    // cout << (vec2) << endl;
  }
  cout << "# [done] (mpz_t, double) sample()" << endl;
  etime = clock();
  secs  = (etime - stime) / (double)CLOCKS_PER_SEC;
  cout << "# [info] (mpz_t, double) took time " << secs << " s" << endl;
#endif

#if 1
  /* mpz_t matrix, mpfr */
  ZZ_mat<mpz_t> M3(dim, dim);
  M3.gen_GM(20 * dim);
  lllReduction(M3, 0.99, 0.51, LM_WRAPPER);
  Z_NR<mpz_t> x3;
  KleinSampler<mpz_t, FP_NR<mpfr_t>> S3(M3);
  FP_NR<mpfr_t> c2, t2;
  c2 = 0.0;
  t2 = 1000.0;

  /* sample Z */
  for (int i = 0; i < 1000; i++)
  {
    x3 = S3.sample_z(c2, t2);
  }
  cout << "# [done] (mpz_t, mpfr_t) sample_z_dgs()" << endl;

  /* sample vec */
  NumVect<Z_NR<mpz_t>> vec3;
  stime = clock();
  for (int i = 0; i < 10000; i++)
  {
    vec3 = S3.sample();
  }
  cout << "# [done] (mpz_t, mpfr_t) sample()" << endl;
  etime = clock();
  secs  = (etime - stime) / (double)CLOCKS_PER_SEC;
  cout << "# [info] (mpz_t, mpfr_t) took time " << secs << " s" << endl;
#endif

  return 0;
}
