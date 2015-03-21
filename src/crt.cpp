/**
*	This program provides an interface for solving the following system of congruence equations:
*	C = c_{i} (mod p_{i})     i = 1,...,k
*	by recombining the c_{i}'s using Chinese Remainder Theorem. 
*
* Compile using the command: g++ crt.cpp -std=c++11 -o crt -lgmp
*
*	Usage:
*	Sample Input: 
*	3			  Number of elements in the array
*	1 3 2		The c_{i}'s or numbers
*	2 5 7		The moduli
*
*	Output:
*	23
*/

#include <iostream>
#include <vector>
#include <gmp.h>
using namespace std;

/**
* Define P = p_1 * p_2 * ... * p_k
* Using Extended Euclidean Algorithm we can find integers r_i, s_i such that
* r_i * p_i + (s_i * P)/p_i = 1 since for each i, p_i and P/p_i are coprime. 
* By substituting e_i = (s_i * P)/p_i, we get one solution to the system of 
* congruence relations as C = c_0 * e_0 + c_1 * e_1 + ... + c_k * e_k
*/
template<class T>
T Recombine_CRT(const T numbers, const T moduli, const mpz_t array_size)
{
	mpz_t *result = new mpz_t[1];
	mpz_t *s = new mpz_t[mpz_get_ui(array_size)];
	mpz_t *t = new mpz_t[mpz_get_ui(array_size)];

	mpz_t gcd;
	mpz_init(gcd);
	mpz_t temp;
	mpz_init(temp);
	mpz_t moduli_lcm;
  	mpz_init(moduli_lcm);
  	
  	mpz_set_ui(moduli_lcm, 1);
  	for (int i = 0; i < mpz_get_ui(array_size); i++)	{
  		mpz_lcm(moduli_lcm, moduli_lcm, moduli[i]);
  	}

  	//Chinese Remainder Algorithm
	for (int i = 0; i < mpz_get_ui(array_size); i++)	{
		mpz_divexact(temp, moduli_lcm, moduli[i]);
  		mpz_gcdext(gcd, s[i], t[i], moduli[i], temp);
  		mpz_mul(s[i], t[i], temp);
  		mpz_addmul(result[0], numbers[i], s[i]);
	}

	mpz_cdiv_r(result[0], result[0], moduli_lcm);

	// Return a non-negative solution
	if (mpz_cmp_ui(result[0], 0) < 0)	{
		mpz_add(result[0], result[0], moduli_lcm);
	}
	
	return result;
}

int main()
{
	//Initializing variables/objects
	mpz_t i;
	mpz_init(i);
	mpz_t array_size;
	mpz_init (array_size);
	mpz_inp_str(array_size, NULL, 10);

	mpz_t *numbers = new mpz_t[mpz_get_ui(array_size)];
	mpz_t *moduli = new mpz_t[mpz_get_ui(array_size)];
	mpz_t *result = new mpz_t[1];

	//Scanning for numbers and moduli arrays as input;
	for (int j = 0; j < mpz_get_ui(array_size); j++)	{
		mpz_inp_str(numbers[j], NULL, 10);
	}
	for (int j = 0; j < mpz_get_ui(array_size); j++)	{
		mpz_inp_str(moduli[j], NULL, 10);
	}
	
  	//Calling the Recombine function
  	result = Recombine_CRT(numbers, moduli, array_size);

  	//Print answer
	mpz_out_str(NULL, 10, result[0]);
	cout << endl;

  	return 0;
}