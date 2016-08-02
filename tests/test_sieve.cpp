#include <cstring>
#include <fplll.h>
#include <../fplll/sieve/sieve_main.h> /* standalone bin */

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

using namespace std;
using namespace fplll;

/**
   @brief Read matrix from `input_filename`.

   @param A
   @param input_filename
   @return
*/
template <class ZT>
void read_matrix(ZZ_mat<ZT> &A, const char *input_filename)
{
  istream *is = new ifstream(input_filename);
  *is >> A;
  delete is;
}


/**
   @brief Read vector from `input_filename` into `b`.

   @param b                vector
   @param input_filename   filename
   @return
*/
template <class ZT>
void read_vector(vector<Z_NR<ZT>> &b, const char *input_filename)
{
  istream *is = new ifstream(input_filename);
  *is >> b;
  delete is;
}


/**
   @brief Test sieve by checking if function returns correct vector.

   @param A              input lattice
   @param b              shortest vector
   @return
*/
template <class ZT>
int test_sieve_alg (ZZ_mat<ZT> &A, IntVect &b, int alg) {
  GaussSieve<ZT, FP_NR<double> > gsieve(A, alg, 0, 0);
  Z_NR<ZT> goal_norm;
  goal_norm = 0;
  gsieve.set_goal_norm2 (goal_norm);
  if (gsieve.alg == 3)
    gsieve.run_3sieve();
  else if (gsieve.alg == 4)
    gsieve.run_4sieve();
  else
    gsieve.run_2sieve();
  NumVect<Z_NR<ZT> > v = gsieve.return_first();
  Z_NR<ZT> tmp;
  Z_NR<ZT> norm_s;
  Z_NR<ZT> norm_b;
  for (int i = 0; i < A.get_cols(); i++)
  {
    tmp.mul(v[i], v[i]);
    norm_s.add(norm_s, tmp);
    tmp.mul(b[i], b[i]);
    norm_b.add(norm_b, tmp);
  }
  if (norm_s != norm_b)
    return 1;
  return 0;
}


template <class ZT>
int test_sieve(ZZ_mat<ZT> &A, IntVect &b) {
  int r = 0;
  r |= test_sieve_alg<ZT>(A, b, 2);
  r |= test_sieve_alg<ZT>(A, b, 3);
  r |= test_sieve_alg<ZT>(A, b, 4);
  return r;
}


/**
   @brief Test sieve for matrix stored in file pointed to by `input_filename`.

   @param input_filename   a path
   @return zero on success
*/
template <class ZT>
int test_filename (const char *input_filename, const char *output_filename) {
  ZZ_mat<ZT> A;
  read_matrix(A, input_filename);
  IntVect b;
  read_vector(b, output_filename);
  return test_sieve<ZT>(A, b);
}

/* 
   Note make check uses the following relative path for the filename.
*/
int main(int /*argc*/, char ** /*argv*/) {

  int status = 0;
  status |= test_filename<mpz_t>(TESTDATADIR "/tests/lattices/example_svp_in",
                                 TESTDATADIR "/tests/lattices/example_svp_out");
  if (status == 0)
  {
    cerr << "All tests passed." << endl;
    return 0;
  }
  else
  {
    return -1;
  }

  return 0;
}
