/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.

This file is part of the fplll Library.

The fplll Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The fplll Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the fplll Library; see the file COPYING.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include "util.h"

using namespace fplll;

template <class ZT>
int llldiff (ZZ_mat<ZT>& B1, ZZ_mat<ZT>& B2, int c, int r)
{
  int test=1, i, j, sg;
  Z_NR<ZT> tmp1;
  Z_NR<ZT> tmp2;

  for (i=0; i<r; i++){
    sg = 1;
    tmp1.abs(B1(i,0));
    tmp2.abs(B2(i,0));
    if (tmp1.cmp(tmp2)!=0){
      //cerr << r << ", 0\n";
      test = 0;
    }
    if (tmp1.cmp(B1(i,0))!=0) sg *=-1;
    if (tmp1.cmp(B2(i,0))!=0) sg *=-1;

    if (sg == 1){
      for (j=1; j<c; j++){
        if (B1(i,j).cmp(B2(i,j))!=0){
          //cerr << i << " " << j << "\n";
          test = 0;
        }
      }
    }
    else{
      for (j=1; j<c; j++){
        tmp1.mul_si(B1(i,j),-1);
        if (tmp1.cmp(B2(i,j))!=0){
          //cerr << i << " " << j << "\n";
          test = 0;
        }
      }
    }
  }

  return (test);
}



/* ********************** */
/*  MAIN **************** */
/* ********************** */

int
main (int argc, char ** argv)
{
  int c, r, ac;
  ZZ_mat<mpz_t> mat1, mat2;

  for (ac = 1; ac < argc; ac++) {
    if (strcmp(argv[ac], "-c") == 0 || strcmp(argv[ac], "-r") == 0) {
      if (ac < argc - 1) ++ac;
    }
    else if (strcmp(argv[ac], "--help") == 0) {
      cout << "Usage: cat matrix1 matrix2 | " << argv[0] << endl;
      return 0;
    }
    else if (argv[ac][0] == '-') {
      cerr << "llldiff: invalid option '" << argv[ac] << "'" << endl;
      return 1;
    }
    else {
      break;
    }
  }

  istream* inputStream;
  if (argv[ac])
    inputStream = new ifstream(argv[ac]);
  else
    inputStream = &cin;

  *inputStream >> mat1;
  *inputStream >> mat2;

  if (argv[ac])
    delete inputStream;

  r = mat1.getRows();
  c = mat1.getCols();

  int difference = !llldiff<mpz_t>(mat1, mat2, c, r);

  if (difference) {
    cerr << "===INVALID RESULT===" << endl;
  }

  return difference;
}
