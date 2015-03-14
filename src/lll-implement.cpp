#include <iostream>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
using namespace std;

//Calculate dot product of a[i] and b[j], used by GSO()
float dotProduct(float a[10][10], int i, float b[10][10], int j, int ineachbasis)
{
    float dot= 0;
    for(int k = 1; k <= ineachbasis; k++)
    {
        dot+=(a[i][k]*b[j][k]);
    }
    return dot;
}

//Apply Gram-Schmidt Process, matrix b gives the Gram-Schmidt basis
void GSO(float a[10][10], float b[10][10], int noofbasis, int ineachbasis, float B[10], float u[10][10])
{
    int i = 1;
    for(int j = 1; j <= ineachbasis; j++) b[i][j] = a[i][j];
    
    B[1] = dotProduct(b, 1, b, 1, ineachbasis);
    
    
    for(i = 2; i<=noofbasis; i++)
    {
        for(int j = 1; j <= ineachbasis; j++) b[i][j] = a[i][j];
        
        for(int j = 1; j <= i-1; j++)
        {
	    //Calculate projection and find orthogonal basis
            u[i][j] = dotProduct(a, i, b, j, ineachbasis)/B[j];
            for(int k = 1; k <= ineachbasis; k++) b[i][k] = b[i][k] - (u[i][j]*b[j][k]);
        }
        
        B[i] = dotProduct(b, i, b, i, ineachbasis);
        
    }
}

//Size Reduction
void RED(float a[10][10], float u[10][10], int k, int l, int noofbasis)
{
        float r = 0.5 + u[k][l];
        int f = floor(r);
        for(int p = 1; p <= noofbasis; p++)
        {
            a[k][p] = a[k][p] - (f*a[l][p]);
        }
        for(int j = 1; j<=l-1; j++)
        {
            u[k][j] = u[k][j] - (f*u[l][j]);
        }
        u[k][l] = u[k][l] - f;
}

//Outputs the result of applying LLL Algorithm on matrix A
void LLLImplement(const Matrix<float> &A)
{
    float a[10][10];
    float b[10][10];
    float u[10][10];
    float B[10];
    int noofbasis, ineachbasis;
    
    noofbasis = A.getRows();
    ineachbasis = A.getCols();
    
    for(int i = 1; i <= noofbasis; i++)
    {
        for(int j = 1; j <= ineachbasis; j++) 
        {
            a[i][j] = A[i][j];
        }
    }
    
    GSO(a, b, noofbasis, ineachbasis, B, u);
    
    int k = 2;
    while(k<=noofbasis)
    {
        for(int j = k-1; j>=1; j--)
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
            for(int p = 1; p <= ineachbasis; p++) swap(a[k][p], a[k-1][p]);
            GSO(a, b, noofbasis, ineachbasis, B, u);
            k = max(2, k-1);
        }
    }
    
    //Output the result
    for(int i = 1; i <= noofbasis; i++)
    {
        for(int j = 1; j <= ineachbasis; j++)
        {
            cout<<a[i][j]<<" ";
        }
        cout<<"\n";
    }
}