/* ---------------------------
   Random generator for mpz_t
   --------------------------- */

#include <time.h>

#ifndef FPLLL_NR_RAND_H
#define FPLLL_NR_RAND_H

FPLLL_BEGIN_NAMESPACE

class RandGen {
public:
  static void init() {
    initialized = true;
    gmp_randinit_default(gmpState);
  }
  static void initWithSeed(unsigned long seed) {
    init();
    gmp_randseed_ui(gmpState, seed);
  }
  static void initWithTime() {
    init();
    gmp_randseed_ui(gmpState, time(NULL));
  }
  static void initWithTime2() {
    init();
    struct timeval time; 
    gettimeofday(&time,NULL);
    gmp_randseed_ui(gmpState, (time.tv_sec * 1000)
                    + (time.tv_usec / 1000));
  }
  static bool getInitialized() {
    return initialized;
  }
  static gmp_randstate_t& getGMPState() {
    if (!initialized) init();
    return gmpState;
  }
  static gmp_randstate_t gmpState;
private:
  static bool initialized;
};


class RandGenInt {
public:
  static void init() {
    initialized = true;
    srand(time(NULL));
  }
  static int get() {
    if (!initialized)
      init();
    return rand();
  }
  static int getBit() {
    if (!initialized)
      init();
    if (rand() % 2 == 0)
      return -1;
    else
      return 1;
  }
  static bool initialized;
};


FPLLL_END_NAMESPACE

#endif
