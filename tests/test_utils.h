#ifndef TEST_UTILS_H
#define TEST_UTILS_H 

#include <fplll.h>
using namespace std;
using namespace fplll;

#define read_file(X, input_filename) {\
  ifstream is;\
  is.exceptions(std::ifstream::failbit | std::ifstream::badbit);\
  try {\
    is.open(input_filename);\
    is >> X;\
    is.close();\
  }\
  catch (const ifstream::failure&) {\
    status = 1;\
    cerr << "Error by reading " << input_filename << "." << endl;\
  }\
}

/**
   @brief Read matrix from `input_filename`.

   @param A matrix
   @param input_filename
   @return zero if the file is correctly read, 1 otherwise.
*/
template <class ZT> int read_matrix(ZZ_mat<ZT> &A, const char *input_filename)
{
  int status = 0;
  read_file(A, input_filename);
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
  read_file(b, input_filename);
  return status;
}

#endif /* TEST_UTILS_H */
