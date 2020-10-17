#pragma once
#include <memory>

/**
 * Searches all points in the lattice spanned by the columns of mu with radius <= initial_radius, and returns the 
 * coefficients (w.r.t. mu) of the shortest, nonzero point among them.
 * 
 * Parameters:
 * mu				The lattice basis, expected to be a upper triangular matrix; in the case of arbitrary lattices, use
 *					the Gram-Schmidt coefficient matrix. Expected to be a pointer to a matrix of size levels x levels in
 *					column-major storage, i.e. the columns of this matrix are stored in successive batches of contigous memory
 * 
 * levels			The rank of the lattice to search
 * 
 * initial_radius	The radius within which to search for lattice points
 */
std::unique_ptr<float[]> search_enumeration_cuda(const float *mu, const unsigned int levels, const float initial_radius)