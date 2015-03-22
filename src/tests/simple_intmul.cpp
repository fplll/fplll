#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<fplll.h>
#include<nr.h>

using namespace std;
using namespace fplll;


template<class T> int zeroIdentity()
{
    T a, b, res;
    a = (long) 0;
    b.randb(30);
    res = a;

    a.mul(a, b);
    
    a = (long) 1;
    res = b;
    
    a.mul(a, b);
    return a.cmp(res);
}

template<class T> int commutative()
{
    T a, b, lhs, rhs;
    a.randb(30);
    b.randb(30);

    lhs.mul(a, b);
    rhs.mul(b, a);

    return lhs.cmp(rhs);
}

template<class T> int associative()
{
    T a, b, c, lhs, rhs;
    a.randb(30);
    b.randb(30);
    c.randb(30);

    lhs.mul(a, b);
    lhs.mul(lhs, c);

    rhs.mul(b, c);
    rhs.mul(rhs, a);

    return lhs.cmp(rhs);
}

template<class T> int distributive()
{
    T a, b, c, lhs, rhs, temp;
    a.randb(30);
    b.randb(30);
    c.randb(30);

    lhs.add(b, c);
    lhs.mul(a, lhs);

    temp.mul(a, b);
    rhs.mul(a, c);
    rhs.add(temp, rhs);
    
    return lhs.cmp(rhs);
}

template<class T> void test_mul()
{
    if( zeroIdentity<T>() )
        FPLLL_ABORT("Zero Identity failed.\n");
    if( commutative<T>() )
        FPLLL_ABORT("Commutative Identity failed.\n");
    if( associative<T>() )
        FPLLL_ABORT("Associative Identity failed.\n");
    if( distributive<T>() )
        FPLLL_ABORT("Distributive Identity failed.\n");
}

int main()
{
    RandGen g;
    g.initWithTime();

    test_mul<Z_NR<long> >();
    test_mul<Z_NR<double> >();
    test_mul<Z_NR<mpz_t> >();
    
    return 0;
}
