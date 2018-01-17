/* ---------------------------
   Random generator for mpz_t
   --------------------------- */

#include <time.h>

#ifndef FPLLL_NR_RAND_H
#define FPLLL_NR_RAND_H

FPLLL_BEGIN_NAMESPACE

class RandGen {
public:
 /**    
	@brief Create a random generator (using gmp):

	@param init             Initialize State for a default algorithm (Mersenne Twister) 
	@param init_with_seed   Set an initial seed value into state. The seed is provided by the user  
        @param init_with_time   Set an initial seed value into state taken from the clock
        @param init_with_time2  As previous, but it applies a function to the time provided by the clock
	@param get_initialized  If an initial seed was set into state, returns True, otherwise it returns False 
	@param get_gmp_state    Returns the current state of the generator 
 */
    static void init() {
    initialized = true;
    gmp_randinit_default(gmp_state);
  }
  static void init_with_seed(unsigned long seed) {
    init();
    gmp_randseed_ui(gmp_state, seed);
  }
  static void init_with_time() {
    init();
    gmp_randseed_ui(gmp_state, time(NULL));
  }
  static void init_with_time2() {
    init();
    struct timeval time; 
    gettimeofday(&time,NULL);
    gmp_randseed_ui(gmp_state, (time.tv_sec * 1000)
                    + (time.tv_usec / 1000));
  }
  static bool get_initialized() {
    return initialized;
  }
  static gmp_randstate_t& get_gmp_state() {
    if (!initialized) init();
    return gmp_state;
  }
  static gmp_randstate_t gmp_state;
private:
  static bool initialized;
};


class RandGenInt {
public:
 /**    
	@brief Create a Random Generator

        @param init  Initialize a random number generator. It does not return somethinng     
        @param get   Returns a pseudo-random integral number in the range between 0 and RAND_MAX (RAND_MAX is at least  32767)
        @get_bit   Returns a random bit 
 */
  static void init() {
    initialized = true;
    srand(time(NULL));
  }
  static int get() {
    if (!initialized)
      init();
    return rand();
  }
  static int get_bit() {
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
