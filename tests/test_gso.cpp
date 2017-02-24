/* Copyright (C) 2015 Martin Albrecht

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
#include <random>
#include <gso.h>
#include <gso_gram.h>
#include <gso_base.h>
#include <nr/matrix.h>


using namespace std;
using namespace fplll;


template <class ZT, class FT> 
bool rs_are_equal(MatGSO<ZT,FT> M1, MatGSOGram<ZT,FT> M2, FT error)
{
  int r = M1.r.get_rows();
  FT maxdiff = 0.0;
  FT diff;
  FT relativation;
  Matrix<FT> r1 = M1.get_r_matrix();
  Matrix<FT> r2 = M2.get_r_matrix();
  for( int i = 0; i < r; i++)
  {
    for( int j = 0; j <= i; j++)
    {
      diff = abs(r1(i,j) - r2(i,j));
      relativation = abs(r1(i,j) + r2(i,j));
      if (!r1(i,j).is_nan() && !r2(i,j).is_nan()) {
		if (relativation > error) {
      		diff = diff/relativation;
      	}
      	if ((diff > error) || (diff.is_nan())) { 
      		cerr << i << " " << j << " " << r1(i,j) << " " << r2(i,j) << "\n"; 
      		r1.print(cerr); 
      		cerr << endl;
      		r2.print(cerr);
      		cerr << endl;
      		return false; 
      	}	
      } else if (r1(i,j).is_nan() && r2(i,j).is_nan()) {

      } else {
      	    cerr << i << " " << j << " " << r1(i,j) << " " << r2(i,j) << "\n"; 
      		r1.print(cerr); 
      		cerr << endl;
      		r2.print(cerr);
      		cerr << endl;
      		return false;
      }

    }
  }
  cerr << "Maximum relative difference: " << maxdiff << endl;
  return true;
}



template <class ZT> void read_matrix(ZZ_mat<ZT> &A, const char *input_filename)
{
  istream *is = new ifstream(input_filename);
  *is >> A;
  delete is;
}

/**
   @brief Test the tester.

   @param A
   @return zero on success.
*/

template <class ZT, class FT> int test_test(ZZ_mat<ZT> &A)
{
  // TEST A
  // Method:
  // - Apply 'normal' MatGSO to A.
  // - Extract r-matrix of A
  // - Compute G = A^T A.
  // - Apply gram MatGSO to G.
  // - Extract r-matrix of G
  // -> The r-matrices should be equal.

  // TEST B
  // Apply some 'random' elementary operation on A and on G
  // (of course, the operations on A and G are only 'abstractly' the same)
  // check if their r-matrices are still equal.

  // TEST A
  // ----------------------------
  int r = A.r;	

  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  MatGSO<Z_NR<ZT>, FP_NR<FT>> Mbuf(A, U, UT, 1);
  Mbuf.update_gso();
  Matrix<Z_NR<ZT>> G = Mbuf.get_g_matrix();

  MatGSO<Z_NR<ZT>, FP_NR<FT>> M(A, U, UT, 0);
  M.update_gso();
  MatGSOGram<Z_NR<ZT>, FP_NR<FT>> M2(G, U, UT, 1);
  M2.update_gso();

  FP_NR<FT> err = .001;
  bool retvalue1 =  rs_are_equal(M, M2, err);

  // TEST B
  // ------------------------
  
  for(int i = 0; i < rand() % 10 + 1; i++) {
  	int k = rand() % r;
  	int j = rand() % r;
  	M.move_row(k,j);
  	M2.move_row(k,j);
  }
  M.update_gso();
  M2.update_gso();
  bool retvalue2 =  rs_are_equal(M, M2, err);

  for(int i = 0; i < rand() % 10 + 1; i++) {
  	int k = rand() % r;
  	int j = rand() % r;
  	M.row_add(k,j);
  	M2.row_add(k,j);
  }
  M.update_gso();
  M2.update_gso();
  bool retvalue3 =  rs_are_equal(M, M2, err);



  return  (!retvalue1)*1 + (!retvalue2)*2 + (!retvalue3)*4;

  //if (retvalue1) { return 1; } 
  //return 0;
}


template <class ZT, class FT>
int test_ggso(ZZ_mat<ZT> &A)
{

  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  //int status = 0;

  // zero on success
  return test_test<ZT,FT>(A);

}


template <class ZT, class FT>
int test_filename(const char *input_filename)
{
  ZZ_mat<ZT> A;
  read_matrix(A, input_filename);
  int retvalue = test_ggso<ZT,FT>(A);
  if (retvalue & 1) 
  {
  	cerr << input_filename << " shows different GSO-outputs for grammatrix representation and basis representation.\n";
  }
  if (retvalue & 2) 
  {
  	cerr << input_filename << " shows different GSO-outputs for grammatrix representation and basis representation after moving rows.\n";
  }
  if (retvalue & 4) 
  {
  	cerr << input_filename << " shows different GSO-outputs for grammatrix representation and basis representation after adding rows.\n";
  }
  if (retvalue > 0) { return 1; } else { return 0; }
}

/**
   @brief Construct d × (d+1) integer relations matrix with bit size b and test LLL.

   @param d                dimension
   @param b                bit size
   @param method           LLL method to test
   @param float_type       floating point type to test
   @param flags            flags to use
   @param prec             precision used for is_lll_reduced

   @return zero on success
*/

template <class ZT, class FT>
int test_int_rel(int d, int b, FloatType float_type = FT_DEFAULT, int prec = 0)
{
  ZZ_mat<ZT> A;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  int retvalue = test_ggso<ZT,FT>(A);
  if (retvalue >= 1)
  {
  	cerr << "Integer relation matrix with parameters " << d << " and " << b << " shows different GSO-outputs for grammatrix representation and basis representation.\n";
  	return 1;
  }
  return 0;
}


int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;

  status |= test_filename<mpz_t,double>("lattices/example2_in");
  status |= test_filename<mpz_t,double>("lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t,double>("lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t,double>("lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t,double>("lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t,double>("lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t,double>(50, 20);
  status |= test_int_rel<mpz_t,double>(40, 10);


  status |= test_filename<mpz_t,mpfr_t>("lattices/example2_in");
  status |= test_filename<mpz_t,mpfr_t>("lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t,mpfr_t>("lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t,mpfr_t>("lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t,mpfr_t>("lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t,mpfr_t>("lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t,mpfr_t>(50, 20);
  status |= test_int_rel<mpz_t,mpfr_t>(40, 10);

#ifdef FPLLL_WITH_LONG_DOUBLE
  status |= test_filename<mpz_t,long double>("lattices/example2_in");
  status |= test_filename<mpz_t,long double>("lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t,long double>("lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t,long double>("lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t,long double>("lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t,long double>("lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t,long double>(50, 20);
  status |= test_int_rel<mpz_t,long double>(40, 10);
#endif
#ifdef FPLLL_WITH_QD
  status |= test_filename<mpz_t,dd_real>("lattices/example2_in");
  status |= test_filename<mpz_t,dd_real>("lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t,dd_real>("lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t,dd_real>("lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t,dd_real>("lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t,dd_real>("lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t,dd_real>(50, 20);
  status |= test_int_rel<mpz_t,dd_real>(40, 10);

  status |= test_filename<mpz_t,qd_real>("lattices/example2_in");
  status |= test_filename<mpz_t,qd_real>("lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t,qd_real>("lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t,qd_real>("lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t,qd_real>("lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t,qd_real>("lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t,qd_real>(50, 20);
  status |= test_int_rel<mpz_t,qd_real>(40, 10);
#endif
#ifdef FPLLL_WITH_DPE
  status |= test_filename<mpz_t,dpe_t>("lattices/example2_in");
  status |= test_filename<mpz_t,dpe_t>("lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t,dpe_t>("lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t,dpe_t>("lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t,dpe_t>("lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t,dpe_t>("lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t,dpe_t>(50, 20);
  status |= test_int_rel<mpz_t,dpe_t>(40, 10);
#endif
  
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
