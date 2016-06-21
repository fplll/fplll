//gcc -std=gnu99 test_gauss_z.c dgs*.c -o test_gauss_z -lm -lmpfr

#include <math.h>
#include "dgs_gauss.h"

#define NTRIALS 1<<18
#define TOLERANCE 0.1
#define BOUND 2

int test_defaults_dp() {
  dgs_disc_gauss_dp_t *self;
  self = dgs_disc_gauss_dp_init(3.0, 0, 6, DGS_DISC_GAUSS_DEFAULT);
  if (self->algorithm != DGS_DISC_GAUSS_UNIFORM_TABLE)
    dgs_die("automatic choice of uniform table algorithm failed (%d)", self->algorithm);
  dgs_disc_gauss_dp_clear(self);

  self = dgs_disc_gauss_dp_init(3.0, 0, 1<<10, DGS_DISC_GAUSS_DEFAULT);
  if (self->algorithm != DGS_DISC_GAUSS_UNIFORM_TABLE)
    dgs_die("automatic choice of uniform table algorithm failed (%d)", self->algorithm);
  dgs_disc_gauss_dp_clear(self);

  self = dgs_disc_gauss_dp_init(3.0, 0, 1<<14, DGS_DISC_GAUSS_DEFAULT);
  if (self->algorithm != DGS_DISC_GAUSS_UNIFORM_LOGTABLE)
    dgs_die("automatic choice of uniform table algorithm failed (%d)", self->algorithm);
  dgs_disc_gauss_dp_clear(self);


  double sigma2 = sqrt(1.0/(2*log(2.0)));
  self = dgs_disc_gauss_dp_init(1024*sigma2, 0, 6, DGS_DISC_GAUSS_DEFAULT);
  if (self->algorithm != DGS_DISC_GAUSS_SIGMA2_LOGTABLE)
    dgs_die("automatic choice of uniform table algorithm failed (%d)", self->algorithm);
  dgs_disc_gauss_dp_clear(self);


  printf("passed\n");
  return 0;
}

int test_uniform_boundaries_dp(double sigma, double c, size_t tau, dgs_disc_gauss_alg_t algorithm) {
  dgs_disc_gauss_dp_t *self = dgs_disc_gauss_dp_init(sigma, c, tau, algorithm);

    printf("σ: %6.2f, c: %6.2f. τ: %2ld, precision: double, algorithm: %d\n", self->sigma, self->c, self->tau, self->algorithm);


  long lower_bound = ((long)self->c) - ceil(self->sigma * self->tau);
  long upper_bound = ((long)self->c) + ceil(self->sigma * self->tau);

  for(size_t i=0; i<NTRIALS; i++) {
    long r = self->call(self);
    if(__DGS_UNLIKELY(r < lower_bound))
      dgs_die("r (%ld) < lower_bound (%ld)", r, lower_bound);
    else if(__DGS_UNLIKELY(r > upper_bound))
      dgs_die("r (%ld) > upper_bound (%ld)", r, upper_bound);
  }
  return 0;
}

/**
 We test if the proportional probabilities holds
*/
#define RHO(x) exp(-(x*x)/(2*sigma*sigma))

int test_ratios_dp(double sigma, size_t tau, dgs_disc_gauss_alg_t algorithm) {
  printf("σ: %6.2f, c:    0.0. τ: %2ld, precision: double, algorithm: %d\n",sigma, tau, algorithm);

  dgs_disc_gauss_dp_t *self = dgs_disc_gauss_dp_init(sigma, 0, tau, algorithm);

  double ctr[2*BOUND+1];

  for(size_t i=0; i<NTRIALS; i++) {
    long r = self->call(self);
    if (abs(r) <= BOUND)
      ctr[r+BOUND] += 1;
  }

  for(long i=-BOUND; i<=BOUND; i++) {
    double left  = ctr[BOUND+1]/ctr[BOUND+i];
    double right = RHO(0)/RHO(i);
    if (fabs(log(left/right)) >= 0.4)
      dgs_die("exp(-((-c)²)/(2σ²))/exp(-((%d-c)²)/(2σ²)) = %7.5f != %7.5f (%7.5f)", i, right, left, fabs(log(left/right)));
  }
  return 0;
}

int test_mean_dp(double sigma, double c, size_t tau, dgs_disc_gauss_alg_t algorithm) {

  printf("σ: %6.2f, c: %6.2f. τ: %2ld, precision: double, algorithm: %d\n",sigma, c, tau, algorithm);

  dgs_disc_gauss_dp_t *self = dgs_disc_gauss_dp_init(sigma, c, tau, algorithm);

  double mean = 0.0;

  for(size_t i=0; i<NTRIALS; i++) {
    long r = self->call(self);
    mean += r;

  }

  mean /=NTRIALS;

  if(fabs(mean - c) > TOLERANCE)
    dgs_die("expected mean %6.2f but got %6.2f",c, mean);

  return 0;
}


int main(int argc, char *argv[]) {
  printf("# testing defaults #\n");
  test_defaults_dp();
  printf("\n");

  printf("# testing proportional probabilities #\n");
  test_ratios_dp( 3.0, 6, DGS_DISC_GAUSS_DEFAULT);
  test_ratios_dp( 3.0, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_ratios_dp( 3.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_ratios_dp( 3.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);
  test_ratios_dp( 2.0, 6, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  test_ratios_dp( 2.0, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_ratios_dp( 2.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_ratios_dp( 2.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);
  test_ratios_dp( 4.0, 3, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  test_ratios_dp( 4.0, 3, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_ratios_dp( 4.0, 3, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_ratios_dp( 4.0, 3, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);
  test_ratios_dp(15.4, 3, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  test_ratios_dp(15.4, 3, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_ratios_dp(15.4, 3, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_ratios_dp(15.4, 3, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);
  printf("\n");

  printf("# testing [⌊c⌋-⌈στ⌉,…, ⌊c⌋+⌈στ⌉] boundaries #\n");
  test_uniform_boundaries_dp( 3.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  test_uniform_boundaries_dp(10.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  test_uniform_boundaries_dp( 3.3, 1.0, 1, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  test_uniform_boundaries_dp( 2.0, 1.5, 2, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  printf("\n");

  test_uniform_boundaries_dp( 3.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_uniform_boundaries_dp(10.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_uniform_boundaries_dp( 3.3, 1.0, 1, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_uniform_boundaries_dp( 2.0, 1.5, 2, DGS_DISC_GAUSS_UNIFORM_TABLE);
  printf("\n");

  test_uniform_boundaries_dp( 3.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_uniform_boundaries_dp(10.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_uniform_boundaries_dp( 3.3, 1.0, 1, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_uniform_boundaries_dp( 2.0, 2.0, 2, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  printf("\n");

  printf("# testing c is center #\n");
  test_mean_dp( 3.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  test_mean_dp(10.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  test_mean_dp( 3.3, 1.0, 6, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  test_mean_dp( 2.0, 1.5, 6, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  printf("\n");

  test_mean_dp( 3.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_mean_dp(10.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_mean_dp( 3.3, 1.0, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);
  test_mean_dp( 2.0, 1.5, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);
  printf("\n");

  test_mean_dp( 3.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_mean_dp(10.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_mean_dp( 3.3, 1.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  test_mean_dp( 2.0, 2.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  printf("\n");

  test_mean_dp( 3.0, 0.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);
  test_mean_dp(10.0, 0.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);
  test_mean_dp( 3.3, 1.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);
  test_mean_dp( 2.0, 2.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);

  printf("\n");

  return 0;
}
