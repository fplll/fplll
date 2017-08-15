#ifndef TEST_UTILS_H
#define TEST_UTILS_H 

#include <fplll.h>
using namespace std;
using namespace fplll;

/**
   @brief Read matrix from `input_filename`.

   @param A matrix
   @param input_filename
   @return zero if the file is correctly read, 1 otherwise.
*/
template <class ZT> int read_matrix(ZZ_mat<ZT> &A, const char *input_filename)
{
  int status = 0;
  ifstream is(input_filename);
  if (!is)
  {
    status = 1;
    cerr << "Could not open file " << input_filename << "." << endl;
  }  // throw std::runtime_error("could not open input file");
  is >> A;
  return status;
}

/**
   @brief Read vector from `input_filename` into `b`.

   @param b                vector
   @param input_filename   filename
   @return zero if the file is correctly read, 1 otherwise.
*/

template <class ZT> int read_vector(vector<Z_NR<ZT>> &b, const char *input_filename)
{
  int status = 0;
  ifstream is;
  is.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    is.open(input_filename);
    is >> b;
    is.close();
  }
  catch (ifstream::failure e) {
    status = 1;
    cerr << "Error by reading " << input_filename << "." << endl;
  }
  return status;
}

#endif /* TEST_UTILS_H */
