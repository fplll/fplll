#ifndef TEST_UTILS_H
#define TEST_UTILS_H 

#include <fplll.h>
using namespace std;
using namespace fplll;

/**
   @brief Read T from `input_filename`.

   @param X T (T is usually a ZZ_mat<ZT> or a vector<Z_NR<ZT>>
   @param input_filename
   @return zero if the file is correctly read, 1 otherwise.
*/
template <class T> int read_file_process(T &X, const char *input_filename) {
  int status = 0;
  ifstream is;
  is.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    is.open(input_filename);
    is >> X;
    is.close();
  }
  catch (const ifstream::failure&) {
    status = 1;
    cerr << "Error by reading " << input_filename << "." << endl;
  }

  return status;
}
#endif /* TEST_UTILS_H */
