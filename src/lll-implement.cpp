#include <iostream>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
using namespace std;

//Calculate dot product of a[i] and b[j], used by GSO()
int dotProduct(const Matrix<int>& a, int i,const Matrix<int>& b, int j, int ineachbasis)
{
    int dot= 0;
    for(int k = 0; k < ineachbasis; k++)
    {
        dot+=(a[i][k]*b[j][k]);
    }
    return dot;
}

//Apply Gram-Schmidt Process, matrix b gives the Gram-Schmidt basis
void GSO(Matrix<int>& a, Matrix<int>& b, int noofbasis, int ineachbasis, int B[], Matrix<float>& u)
{
    int i = 0;
    for(int j = 0; j < ineachbasis; j++) b[i][j] = a[i][j];
    
    B[0] = dotProduct(b, 0, b, 0, ineachbasis);
    
    
    for(i = 1; i<noofbasis; i++)
    {
        for(int j = 0; j < ineachbasis; j++) b[i][j] = a[i][j];
        
        for(int j = 0; j <= i-1; j++)
        {
	    //Calculate projection and find orthogonal basis
            u[i][j] = dotProduct(a, i, b, j, ineachbasis)/B[j];
            for(int k = 0; k < ineachbasis; k++) b[i][k] = b[i][k] - (u[i][j]*b[j][k]);
        }
        
        B[i] = dotProduct(b, i, b, i, ineachbasis);
        
    }
}

//Size Reduction
void RED(Matrix<int>& a, Matrix<float>& u, int k, int l, int noofbasis)
{
        float r = 0.5 + u[k][l];
        int f = floor(r);
        for(int p = 0; p < noofbasis; p++)
        {
            a[k][p] = a[k][p] - (f*a[l][p]);
        }
        for(int j = 0; j<l-1; j++)
        {
            u[k][j] = u[k][j] - (f*u[l][j]);
        }
        u[k][l] = u[k][l] - f;
}

//Outputs the result of applying LLL Algorithm on matrix A
void LLLImplement(const Matrix<int> &A)
{
    int noofbasis, ineachbasis;
    
    noofbasis = A.getRows();
    ineachbasis = A.getCols();
    
    Matrix<int> a(noofbasis, ineachbasis);
    Matrix<float> u(noofbasis, ineachbasis);
    Matrix<int> b(noofbasis, ineachbasis);
    int B[noofbasis];
    
    for(int i = 0; i < noofbasis; i++)
    {
        for(int j = 0; j < ineachbasis; j++) 
        {
            a[i][j] = A[i][j];
        }
    }
    
    GSO(a, b, noofbasis, ineachbasis, B, u);
    
    int k = 1;
    while(k<noofbasis)
    {
        for(int j = k-1; j>=0; j--)
        {
	    //Length reduce a[k] and calculate correct u[k][j] values
            if(u[k][j]>0.5 || u[k][j]<-0.5)
                RED(a, u, k, j, noofbasis);
        }
        if(B[k] >= (0.75 - (u[k][k-1]*u[k][k-1]))*B[k-1] )	//Lovasz Condition (delta) = 0.75
        {
            k++;
        }
        else	//Lovasz Condition fails, swap rows and calculate new values
        {
            for(int p = 0; p < ineachbasis; p++) swap(a[k][p], a[k-1][p]);
            GSO(a, b, noofbasis, ineachbasis, B, u);
            k = max(1, k-1);
        }
    }
    
    //Output the result
    for(int i = 0; i < noofbasis; i++)
    {
        for(int j = 0; j < ineachbasis; j++)
        {
            cout<<a[i][j]<<" ";
        }
        cout<<"\n";
    }
    a.clear();
    b.clear();
    u.clear();
}