/**
*	This program provides an interface for solving the following system of congruence equations:
*	C = c_{i} (mod p_{i})     i = 1,...,k
*	by recombining the c_{i}'s using Chinese Remainder Theorem. 
*
* Compile using the command: g++ crt.cpp -std=c++11 -o crt
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

using namespace std;

// Structure "egcd" to store the GCD, and Bezout Coefficients such that a*x + b*y = gcd(a, b)
typedef struct {
  int gcd;
  int x;
  int y;  
} egcd;

// Recursive function for Extended Euclid's Algorithm
template<class T>
egcd Extended_Euclid(T a, T b)
{
  // Swap elements if a is greater than b
  if ( a > b) {
    T t = a;
    a = b;
    b = t;
  }

  egcd e;
  if ( a == 0)  {
      e = {b, 0, 1};  
    return e;
  }
  else  {
    e = Extended_Euclid(b%a, a);
    T temp = e.y;
    e.y = e.x;
    e.x = temp; 
    e = {e.gcd, e.x-((b/a)*e.y), e.y};
    return e;
  }
}

/**
* Define P = p_1 * p_2 * ... * p_k
* Using Extended Euclidean Algorithm we can find integers r_i, s_i such that
* r_i * p_i + (s_i * P)/p_i = 1 since for each i, p_i and P/p_i are coprime. 
* By substituting e_i = (s_i * P)/p_i, we get one solution to the system of 
* congruence relations as C = c_0 * e_0 + c_1 * e_1 + ... + c_k * e_k
*/
template<class T>
T Recombine_CRT(vector<T>& numbers, vector<T>& moduli)
{
  if (numbers.size() == 0 || moduli.size() == 0)  {
    cout << "Arguments to function should be non-empty lists." << endl;
    return -1;
  }
  if (numbers.size() != moduli.size())  {
    cout << "Arguments to function should be lists of same length." << endl;
    return -1;
  }

  T C = 0;
  T P = moduli[0];
  egcd e;
  vector<int> E(numbers.size());
  for ( int i = 1; i< moduli.size(); i++) {
    P = P*moduli[i];
  }
  
  for ( int i = 0; i < moduli.size(); i++) {
    e = Extended_Euclid(moduli[i], P/moduli[i]);
    E[i] = (e.y*P)/moduli[i];
    C = C + numbers[i]*E[i]; 
  }

  //To find a non-negative solution
  if( C < 0)  {
    C = P + C;
  }
  return C%P;
}

int main()
{
	//Initializing variables/objects
  int array_size;
  int result;
  cin >> array_size;
  vector<int> numbers(array_size);
  vector<int> moduli(array_size);
     
	//Scanning for c_{i} and p_{i} as input;
  for (int i=0;i<array_size;i++) {
    	cin >> numbers[i];
	}
	for (int i=0;i<array_size;i++) {
	   	cin >> moduli[i];
	}

	//Calling the Recombine function
  result = Recombine_CRT(numbers, moduli);

  //Print answer
	cout << result << endl;	

  return 0;
}