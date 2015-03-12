#include "matrix.h"
#include <math.h>

// The function takes a non-orthogonal matrix m, which is orthogonalised using Gram-Schimdt method

void gramSchmidt(const Matrix<float>& m, Matrix<float>& n)	
{   
    float dotprod;
    int nbasis = m.getRows();
    int dim = m.getCols();
    
    for(int i=0; i<nbasis; ++i)
    { 
        for(int j=0; j<dim; j++)
        {
        //Initialize the vector
            n[i][j] = m[i][j]; 
        }
        for(int k=i+1; k<nbasis; k++)
        {
            dotprod = 0;
            //compute the dot product
            for(int j=0; j<dim; j++)
            {
                dotprod+=(m[k][j]*n[i][j]);
            }
            //Substract the product of dot product with the current set of orthogonal vectors
            for(int j=0; j<dim; j++)
            {
                m[k][j]-=(dotprod*m[i][j]);
            }
        }
    }
}
