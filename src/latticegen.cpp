/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.


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

#include "config.h"
#include "util.h"

using namespace fplll;

void printHelp() {
  cout << "Usage: latticegen [-randseed [<int> | 'time']] options\n"
       << "Options : "<< endl
       << " r <d> <b> : gen_intrel"<<endl
       << " s <d> <b> <b2> : gen_simdioph"<< endl
       << " u <d> <b> : gen_uniform"<< endl
       << " n <d> <b> <q> : gen_ntrulike"<<endl
       << " N <d> <b> <q> : gen_ntrulike2"<<endl
       << " a <d> <f> : gen_ajtai"<<endl
       << " A <d> : gen_ajtai2"<<endl;
}

void printVersion() {
  cout << "latticegen (fplll) " << VERSION << endl
        << "Copyright 2005-2011 Damien Stehle, David Cade, Xavier Pujol." << endl
        << "fplll is free software. You can redistribute it and/or modify" << endl
        << "it under the terms of the GNU Lesser General Public License as published by" << endl
        << "the Free Software Foundation, either version 2.1 of the License, or" << endl
        << "(at your option) any later version." << endl;
}

void fatalError(const char* message) {
  cerr << "latticegen: " << message << endl
       << "Try 'latticegen --help' for more information" << endl;
  exit(1);
}

int main(int argc, char* argv[]) {
  int iArg = 1;
  
  if (argc - iArg < 1 || strcmp(argv[iArg], "--help") == 0) {
    printHelp();
    return 0;
  }
  else if (strcmp(argv[iArg], "--version") == 0) {
    printVersion();
    return 0;
  }
  else if (strcmp(argv[iArg], "-randseed") == 0) {
    iArg++;
    if (argc - iArg < 1) fatalError("option '-randseed' requires an argument");
    if (strcmp(argv[iArg], "time") == 0) {
      RandGen::initWithTime();
    }
    else {
      RandGen::initWithSeed(atol(argv[iArg]));
    }
    iArg++;
  }
  
  if (argc - iArg < 2) fatalError("you must specify a method and a dimension");
  char genMethod = argv[iArg][0];
  int d = atoi(argv[iArg + 1]);
  iArg += 2;
  
  IntMatrix m;
  
  //initialization+filling of the matrix
  switch (genMethod)
    {
    case 'r':
      {
	if (argc - iArg < 1) fatalError("method 'r' requires 2 arguments");
	int b=atoi(argv[iArg]);
	m.resize(d, d + 1);
	m.gen_intrel(b);
	break;
      }
    case 's':
      {
	if (argc - iArg < 2) fatalError("method 's' requires 3 arguments");
	int b=atoi(argv[iArg]);
	int b2=atoi(argv[iArg + 1]);
	m.resize(d + 1, d + 1);
	m.gen_simdioph(b, b2);
	break;
      }
    case 'u':
      {
	if (argc - iArg < 1) fatalError("method 'u' requires 2 arguments");
	int b=atoi(argv[iArg]);
	m.resize(d, d);
	m.gen_uniform(b);
	break;
      }
    case 'n':
      {
	if (argc - iArg < 2) fatalError("method 'n' requires 3 arguments");
	int b=atoi(argv[iArg]);
	int q=atoi(argv[iArg + 1]);
	m.resize(2 * d, 2 * d);
	m.gen_ntrulike(b, q);
	break;
      }
    case 'N':
      {
	if (argc - iArg < 2) fatalError("method 'N' requires 3 arguments");
	int b=atoi(argv[iArg]);
	int q=atoi(argv[iArg + 1]);
        m.resize(2 * d, 2 * d);
	m.gen_ntrulike2(b, q);
	break;
      }
    case 'a':
      {
	if (argc - iArg < 1) fatalError("method 'a' requires 2 arguments");
	double alpha=atof(argv[iArg]);
        m.resize(d, d);
	m.gen_ajtai(alpha);
	break;
      }
    case 'A':
      {
	FP_NR<mpfr_t>* w = new FP_NR<mpfr_t>[d];

	for (int i=0; i<d; i++)
	  mpfr_inp_str(w[i].getData(), stdin, 10, GMP_RNDN);

	m.resize(d, d);
	m.gen_ajtai2(w);

	delete[] w;

	break;
      }
    default:
      {
	fatalError("invalid method");
	break;
      }
    }
  cout << m << endl;
  return 0;
}
