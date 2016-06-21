/* Copyright (C) 2016 Martin R. Albrecht

   This file is part of fplll. fplll is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   fplll is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with fplll. If not, see <http://www.gnu.org/licenses/>. */

#include <cstring>
#include <fplll.h>

using namespace std;
using namespace fplll;

template<class FT> int test_arithmetic() {
  FT a, b, c;
  a = 6.0;
  b = 3.0;
  c =-2.0;

  FT d = (a*b + 2*c + 2.0)/c;
  return d.cmp(-8);
}

int main(int argc, char *argv[]) {

  int status = 0;
  status |= test_arithmetic<FP_NR<double> >();
  status |= test_arithmetic<FP_NR<mpfr_t> >();

  if (status == 0) {
    cerr << "All tests passed." << endl;
    return 0;
  } else {
    return -1;
  }

  return 0;
}
