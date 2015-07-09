#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<fplll.h>
#include<nr.h>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;
using namespace fplll;

template<class T> int zeroIdentity() {
    T a, b, res;
    a = (long) 0;
    b.randb(20);
    res = a;

    a.mul(a, b);
    
    a = (long) 1;
    res = b;
    
    a.mul(a, b);
    return a.cmp(res);
}

template<class T> int commutative() {
    T a, b, lhs, rhs;
    a.randb(20);
    b.randb(20);

    lhs.mul(a, b);
    rhs.mul(b, a);

    return lhs.cmp(rhs);
}

template<class T> int associative() {
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

template<class T> int distributive() {
    T a, b, c, lhs, rhs, temp;
    a.randb(20);
    b.randb(20);
    c.randb(20);

    lhs.add(b, c);
    lhs.mul(a, lhs);

    temp.mul(a, b);
    rhs.mul(a, c);
    rhs.add(temp, rhs);
    
    return lhs.cmp(rhs);
}

TEST_CASE("Integer Multiplication", "intmul") {
    RandGen g;
    g.initWithTime();

    SECTION("long") { 
        CHECK(zeroIdentity<Z_NR<long> >() == 0);
        CHECK(commutative<Z_NR<long> >() == 0);
        CHECK(associative<Z_NR<long> >() == 0);
        CHECK(distributive<Z_NR<long> >() == 0);
    }

    SECTION("double") {
        CHECK(zeroIdentity<Z_NR<double> >() == 0);
        CHECK(commutative<Z_NR<double> >() == 0);
        CHECK(associative<Z_NR<double> >() == 0);
        CHECK(distributive<Z_NR<double> >() == 0);
    }

    SECTION("mpz_t") {
        CHECK(zeroIdentity<Z_NR<mpz_t> >() == 0);
        CHECK(commutative<Z_NR<mpz_t> >() == 0);
        CHECK(associative<Z_NR<mpz_t> >() == 0);
        CHECK(distributive<Z_NR<mpz_t> >() == 0);
    }
}
