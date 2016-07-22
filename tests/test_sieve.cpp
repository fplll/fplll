#include <cstring>
#include <fplll.h>

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
int test_svp(ZZ_mat<ZT> &A, IntVect &b) {
  
  IntVect solution;
  IntMatrix u;

  int status = shortest_vector(A, sol_coord, SVPM_PROVED, SVP_DEFAULT);

  if (status != RED_SUCCESS)
  {
    cerr << "Failure: " << get_red_status_str(status) << endl;
    return status;
  }

  vector_matrix_product(sol_coord2, sol_coord, u);
  vector_matrix_product(solution, sol_coord, A);

  Z_NR<ZT> tmp;
  Z_NR<ZT> norm_s;
  Z_NR<ZT> norm_b;

  for (int i = 0; i < A.get_cols(); i++)
  {
    tmp.mul(solution[i], solution[i]);
    norm_s.add(norm_s, tmp);

    tmp.mul(b[i], b[i]);
    norm_b.add(norm_b, tmp);
  }
  if (norm_s != norm_b)
    return 1;

  return 0;
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


int main(int /*argc*/, char ** /*argv*/) {

  int status = 0;
  status |= test_filename<mpz_t>("lattices/dim55_in");
  status |= test_filename<mpz_t>("lattices/example_svp_in",
                                 "lattices/example_svp_out");
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
