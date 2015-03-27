#include <iostream>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
using namespace std;

//Calculate dot product of a[i] and b[j] vectors
float dotProduct(const Matrix<int>& a, int i, const Matrix<int>& b, int j)
{
    float dot= 0;
    int ineachbasis = a.getCols();
    for(int k = 0; k < ineachbasis; k++)
    {
        dot+=(a[i][k]*b[j][k]);
    }
    return dot;
}

float dotProduct(const Matrix<int>& a, int i, const Matrix<float>& b, int j)
{
    float dot= 0;
    int ineachbasis = a.getCols();
    for(int k = 0; k < ineachbasis; k++)
    {
        dot+=(a[i][k]*b[j][k]);
    }
    return dot;
}

//Apply Gram-Schmidt Process, matrix b gives the Gram-Schmidt basis
void IncGSO(Matrix<float>& b, Matrix<float>& u, float B[], const Matrix<int> a, const int k, int kmax)
{
    int noofbasis, ineachbasis;
    noofbasis = a.getRows();
    ineachbasis = a.getCols();
  
    if(k<=kmax) return;
    kmax = k;
    for(int j = 0; j < ineachbasis; j++) b[k][j] = a[k][j];
    for(int j=0; j<=k-1; j++)
    {
        u[k][j] = dotProduct(a, k, b, j, ineachbasis)/B[j];
        for(int p = 0; p < ineachbasis; p++) b[k][p] -= (u[k][j]*b[j][p]);
    }
    B[k] = dotProduct(b, k, b, k, ineachbasis);
}

//Size Reduction
void RED(Matrix<int>& a, Matrix<float>& u, Matrix<int>& H, const int k, const int l)
{
    int ineachbasis = a.getCols();
    if(u[k][l]<0.5 && u[k][l]>-0.5) return;	//check condition for size reduction
    float r = 0.5+u[k][l];
    int q = floor(r);
    for(int p = 0; p < ineachbasis; p++) a[k][p] -= (q*a[l][p]);
    for(int p = 0; p < ineachbasis; p++) H[k][p] -= (q*H[l][p]);
    u[k][l] -= q;
    for(int i=0; i<=l-1; i++)
    {
        u[k][i] -= (q*u[l][i]);
    }
}

//When Lovasz condition fails
void SWAP(Matrix<int>& a, Matrix<float>& u, Matrix<int>& H, Matrix<float>& b, float B[], const int k, const int kmax)
{
    int ineachbasis = a.getCols();
    for(int p = 0; p < ineachbasis; p++)	//swap a[k], a[k-1] vectors
    {
        swap(a[k][p], a[k-1][p]);
        swap(H[k][p], H[k-1][p]);
    }
    for(int j=0; j<=k-2; j++)
    {
        swap(u[k][j], u[k-1][j]);
    }
    float tu = u[k][k-1];	//temporary u
    float tB = B[k] + (tu*tu*B[k-1]);	//temporary B
    u[k][k-1] = (tu*B[k-1])/tB;
    float tb[ineachbasis];
    for(int p = 0; p < ineachbasis; p++) tb[p] = b[k-1][p];
    for(int p = 0; p < ineachbasis; p++) b[k-1][p] = b[k][p] + (tu*tb[p]);
    for(int p = 0; p < ineachbasis; p++) b[k][p] = ((B[k]/tB)*tb[p]) - (u[k][k-1]*b[k][p]);
    B[k] = B[k-1]*B[k]/tB;
    B[k-1] = tB;
    for(int i = k+1; i<=kmax; i++)
    {
        int t = u[i][k];
        u[i][k] = u[i][k-1] - (tu*t);
        u[i][k-1] = t + (u[k][k-1]*u[i][k]);
    }
}

void LLLImplement(Matrix<int>& A)
{
    int noofbasis, ineachbasis;
    
    noofbasis = A.getRows();
    ineachbasis = A.getCols();
    
    Matrix<int> a(noofbasis, ineachbasis);
    Matrix<float> u(noofbasis, ineachbasis);
    Matrix<float> b(noofbasis, ineachbasis);
    Matrix<int> H(noofbasis, ineachbasis);
    float B[noofbasis];
    int k, kmax;
    
    for(int i = 0; i < noofbasis; i++)
    {
        for(int j = 0; j < ineachbasis; j++) 
        {
            a[i][j] = A[i][j];
        }
    }
    //initialize
    k = 1;
    kmax = 0;
    for(int p = 0; p < ineachbasis; p++) b[0][p] = a[0][p];
    B[0] = dotProduct(a, 0, a, 0, ineachbasis);
    for(int i = 0; i < noofbasis; i++)
	for(int j = 0; j < ineachbasis; j++)
	    if(i==j)    H[i][j] = 1;
	    else        H[i][j] = 0;
    
    
    while(k<noofbasis)
    {
        IncGSO(b, u, B, a, k, kmax);	//Incremental Gram-Schmidt
        RED(a, u, H, k, k-1);
        float chck = (0.75 - u[k][k-1]*u[k][k-1]);
        if(B[k] < chck*B[k-1])	//Lovasz Condition
        {
            SWAP(a, u, H, b, B, k, kmax);
            k = max(1, k-1);
        }
        else
        {
            for(int l = k-2; l>=0; l--)
            {
                RED(a, u, H, k, l);
            }
            k++;
        }
    }
    cout<<"\nLLL-Reduced Basis:\n";
    for(int i = 0; i < noofbasis; i++)
    {
        for(int j = 0; j < ineachbasis; j++)
            cout<<a[i][j]<<" ";
        cout<<"\n";
    }
    cout<<"\nTransformation Matrix:\n";
    for(int i = 0; i < noofbasis; i++)
    {
        for(int j = 0; j < ineachbasis; j++)
            cout<<H[i][j]<<" ";
        cout<<"\n";
    }
    
    return 0;
}