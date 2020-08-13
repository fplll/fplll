#include <cstring>
#include <fplll/fplll.h>

using namespace fplll;

int test_reset()
{
  WholeTreeCounter wc;
  LevelTreeCounter lc;
  wc.reset();
  lc.reset();
  return (wc != (uint64_t)0) | (lc != (uint64_t)0);
}

int test_validate()
{
  WholeTreeCounter wc;
  LevelTreeCounter lc;
  wc.invalidate();
  lc.invalidate();
  int status = wc.is_valid() | lc.is_valid();
  wc         = 1;
  std::array<uint64_t, FPLLL_MAX_ENUM_DIM> a;
  a[0] = 1;
  lc   = a;
  return status | (!wc.is_valid()) | (!lc.is_valid());
}

int test_get_nodes()
{
  WholeTreeCounter wc;
  LevelTreeCounter lc;

  wc = 1;
  std::array<uint64_t, FPLLL_MAX_ENUM_DIM> a{};
  a[0] = 1;
  lc   = a;

  int status = (wc.get_nodes() != 1) | (lc.get_total_nodes() != 1);
  status |= (wc.get_nodes() == 3) | (lc.get_total_nodes() == 3);
  status |= (wc.get_total_nodes() == 3);
  status |= (lc.get_nodes() != a);
  wc.invalidate();
  lc.invalidate();
  status |= wc.get_total_nodes() != lc.get_total_nodes();
  wc.reset();
  lc.reset();
  status |= wc.get_total_nodes() != lc.get_total_nodes();
  return status;
}

template <class FT> int test_enum(size_t d)
{
  /**
   * This is stolen from test_enum.cpp
   */

  RandGen::init_with_seed(0x1337);
  ZZ_mat<mpz_t> A = ZZ_mat<mpz_t>(100, 100);
  A.gen_qary_withq(50, 7681);
  lll_reduction(A);
  ZZ_mat<mpz_t> U;
  MatGSO<Z_NR<mpz_t>, FP_NR<FT>> M(A, U, U, 0);
  M.update_gso();

  // Here we just duplicate everything to prevent any interactions.
  FastEvaluator<FP_NR<FT>> evaluator;
  FastEvaluator<FP_NR<FT>> evaluator2;
  Enumeration<Z_NR<mpz_t>, FP_NR<FT>, WholeTreeCounter> enum_obj(M, evaluator);
  Enumeration<Z_NR<mpz_t>, FP_NR<FT>, LevelTreeCounter> enum_obj2(M, evaluator2);

  // Now just enumerate over the lattice, but only do so over the first d many entries.
  FP_NR<FT> max_dist;
  M.get_r(max_dist, 0, 0);
  max_dist *= 0.99;
  enum_obj.enumerate(0, d, max_dist, 0);
  enum_obj2.enumerate(0, d, max_dist, 0);
  // Here we assume enumeration works, and we only care that their nodes are the same
  // across two different counters
  return enum_obj.get_total_nodes() != enum_obj2.get_total_nodes();
}

int main(int argc, char *argv[])
{

  int status = 0;
  status |= test_reset();
  status |= test_validate();
  status |= test_get_nodes();
  status |= test_enum<double>(30);
  if (status == 0)
  {
    std::cerr << "All tests passed." << std::endl;
    return 0;
  }
  else
  {
    return -1;
  }
}
