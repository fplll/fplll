#include <fplll/fplll.h>

using namespace fplll;

template <class T> bool test_ceil()
{

  bool status = false;

  // Case 1: works with fixed values.
  FP_NR<T> value(3.5);
  FP_NR<T> ceiled{};
  ceiled.ceil(value);

  status |= (ceiled != 4);

  // Case 2: calling ceil either gives us a number >= value.
  value = rand();
  ceiled.ceil(value);
  status |= (ceiled < value);
  return status;
}

int main(int, char **)
{
  int status = 0;
  status += test_ceil<mpfr_t>();
  status += test_ceil<double>();
#ifdef FPLLL_WITH_LONG_DOUBLE
  status += test_ceil<long double>();
#endif

#ifdef FPLLL_WITH_QD
  status += test_ceil<dd_real>();
  status += test_ceil<qd_real>();
#endif

#ifdef FPLLL_WITH_DPE
  status += test_ceil<dpe_t>();
#endif

  if (status == 0)
  {
    std::cerr << "All tests passed" << std::endl;
    return 0;
  }
  else
  {
    return -1;
  }
  return 0;
}
