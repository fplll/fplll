#include<iostream>
#include<time.h>
#include<fplll.h>
#include<matrix.h>
#include<matrix.cpp>
#include<nr.h>
#include<lll.h>
#include<lll.cpp>
#include<gso.h>
#include<gso.cpp>
#include<wrapper.h>

using namespace std;
using namespace fplll;

void gso(Matrix<Integer>& Z, Matrix<Float>& Res,
                Matrix<Float>& coeff)
{
    int i, j, k;
    int rows = Z.getRows(), columns = Z.getCols();
    Float tmp;

    IntMatrix U, Uinvt;
    MatGSO<Integer, Float> G(Z,U,Uinvt,0);
    
    G.updateGSO();
    coeff = G.getMuMatrix();

    for(i=0;i<rows;i++)
    {
        for(j=0;j<columns;j++)
        {
            Res(i, j).set_z(Z(i, j));
            for(k=0;k<i;k++)
            {
                tmp.mul( coeff(i, k) , Res(k, j) );
                Res(i, j).sub( Res(i, j), tmp );
            }
        }
    }
}

int test_sizeReduced(Matrix<Float>& coeff)
{
    int i, j, fail;
    int rows = coeff.getRows();
    Float half;
    half = (double) 0.5;

    fail = 0;
    for(i=0;i<rows;i++)
        for(j=0;j<i;j++)
            if( coeff(i, j).cmp(half) > 0 )
                fail = 1;

    return fail;
}

int test_LovaszCondition(Matrix<Float>& Res, Matrix<Float>& coeff,
                double delta)
{
    int j, k;
    int rows = Res.getRows(), columns = Res.getCols();

    vector<Float> mag;
    vector<Float> mu;

    Float tmp1, tmp2;
    for(k = 0; k < rows; k++)
    {
        tmp1 = 0;
        for(j = 0; j < columns; j++)
        {
            tmp2.mul(Res[k][j], Res[k][j]);
            tmp1.add(tmp1, tmp2);
        }
        mag.push_back(tmp1);
    }

    Float lhs, rhs;
    int fail = 0;
    tmp2 = delta;
    for(k = 2;k < rows; k++)
    {
        tmp1.mul( coeff(k, k-1), coeff(k, k-1) );
        tmp1.sub( tmp2, tmp1 );
        lhs.mul( tmp1, mag[k-1] );

        rhs = mag[k];
        if( lhs.cmp(rhs) > 0 )
            fail = 1;
    }

    return fail;
}

int main()
{
    int rows, columns;
    double delta = 0.99;
    srand(time(NULL));

    rows = rand()%10 + 1;
    columns = rows + (rand()%(15 - rows + 1));

    IntMatrix M(rows, columns), U, Uinvt;
    Matrix<Float> Res(rows, columns), coeff(rows, columns);
    M.gen_uniform(4);

    MatGSO<Z_NR<IntegerT>, FP_NR<FloatT> > G(M, U, Uinvt, 0);
    LLLReduction<Z_NR<IntegerT>, FP_NR<FloatT> > L(G, delta, 0.51, 0);
    L.lll();

    gso(G.b, Res, coeff);
    
    int var1, var2;
    var1 = test_sizeReduced(coeff);
    var2 = test_LovaszCondition(Res, coeff, delta);

    FPLLL_CHECK( var1 == 0 , "Size-Reduce condition failed for LLL.\n");
    FPLLL_CHECK( var2 == 0 , "Lovasz condition failed for LLL.\n");

    return 0;
}
