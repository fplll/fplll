/****************************
 *  Z=double specialization
 ****************************/

#ifndef FPLLL_NR_Z_D_H
#define FPLLL_NR_Z_D_H

FPLLL_BEGIN_NAMESPACE


/* specialization double as z */
#ifdef FPLLL_WITH_ZDOUBLE

template<>
inline Z_NR<double>::Z_NR() {}

template<>
inline Z_NR<double>::Z_NR(const Z_NR<double>& z) : data(z.data) {}

template<>
inline Z_NR<double>::~Z_NR() {}

/** get data */
template<>
inline double Z_NR<double>::get_d() const {
  return data;
}

#ifdef FPLLL_WITH_LONG_DOUBLE
template<>
inline long double Z_NR<double>::get_ld() const {
  return static_cast<long double>(data);
}
#endif

template<>
inline long Z_NR<double>::get_si() const {
  return static_cast<long>(rint(data));
}

template<>
inline void Z_NR<double>::get_mpz(mpz_t r) const {
  mpz_set_d (r, data);
}

template<>
inline long Z_NR<double>::exponent() const {
  int intExpo;
  frexp(data, &intExpo);
  return static_cast<long>(intExpo);
}

/** set data */
template<>
inline void Z_NR<double>::set_str(const char* s) {
  data = atof(s);
}

/** comparison */
template<>
inline int Z_NR<double>::cmp(const Z_NR<double>& m) const {
  if (data > m.data) return 1;
  if (data == m.data) return 0;
  return -1;
}

template<>
inline int Z_NR<double>::sgn() const {
  if (data > 0.0) return 1;
  if (data == 0.0) return 0;
  return -1;
}

/** operator */
template<>
inline void Z_NR<double>::operator=(const Z_NR<double>& a) {
  data = a.data;
}

template<>
inline void Z_NR<double>::operator=(long a) {
  data = static_cast<double>(a);
}

template<>
inline bool Z_NR<double>::operator<(const Z_NR<double>& a) const {
  return data < a.data;
}

template<>
inline bool Z_NR<double>::operator<(long a) const {
  return data < static_cast<double>(a);
}

template<>
inline bool Z_NR<double>::operator>(const Z_NR<double>& a) const {
  return data > a.data;
}

template<>
inline bool Z_NR<double>::operator>(long a) const {
  return data > static_cast<double>(a);
}

template<>
inline bool Z_NR<double>::operator<=(const Z_NR<double>& a) const {
  return data <= a.data;
}

template<>
inline bool Z_NR<double>::operator<=(long a) const {
  return data <= static_cast<double>(a);
}

template<>
inline bool Z_NR<double>::operator>=(const Z_NR<double>& a) const {
  return data >= a.data;
}

template<>
inline bool Z_NR<double>::operator>=(long a) const {
  return data >= static_cast<double>(a);
}

template<>
inline bool Z_NR<double>::operator==(const Z_NR<double>& a) const {
  return data == a.data;
}

template<>
inline bool Z_NR<double>::operator==(long a) const {
  return data == static_cast<double>(a);
}

template<>
inline bool Z_NR<double>::operator!=(const Z_NR<double>& a) const {
  return data != a.data;
}

template<>
inline bool Z_NR<double>::operator!=(long a) const {
  return data != static_cast<double>(a);
}

/** arithmetic */
template<>
inline void Z_NR<double>::add(const Z_NR<double>& a, const Z_NR<double>& b) {
  data = a.data + b.data;
}

template<>
inline void Z_NR<double>::add_ui(const Z_NR<double>& a, unsigned int b) {
  data = a.data + static_cast<double>(b);
}

template<>
inline void Z_NR<double>::sub(const Z_NR<double>& a, const Z_NR<double>& b) {
  data = a.data - b.data;
}

template<>
inline void Z_NR<double>::sub_ui(const Z_NR<double>& a, unsigned int b) {
  data = a.data - static_cast<double>(b);
}

template<>
inline void Z_NR<double>::neg(const Z_NR<double>& a) {
  data = -a.data;
}

template<>
inline void Z_NR<double>::mul(const Z_NR<double>& a, const Z_NR<double>& b) {
  data = a.data * b.data;
}

template<>
inline void Z_NR<double>::mul_si(const Z_NR<double>& a, long int b) {
  data = a.data * static_cast<double>(b);
}

template<>
inline void Z_NR<double>::mul_ui(const Z_NR<double>& a, long unsigned int b) {
  data = a.data * static_cast<double>(b);
}

template<>
inline void Z_NR<double>::mul_2si(const Z_NR<double>& a, long b) {
  data = ldexp(a.data, b);
}

template<>
inline void Z_NR<double>::div_2si(const Z_NR<double>& a, long b) {
  data = ldexp(a.data, -b);
}

template<>
inline void Z_NR<double>::addmul(const Z_NR<double>& a, const Z_NR<double>& b) {
  data += a.data * b.data;
}

template<>
inline void Z_NR<double>::addmul_ui(const Z_NR<double>& a, unsigned long b) {
  data += a.data * static_cast<double>(b);
}

template<>
inline void Z_NR<double>::addmul_si(const Z_NR<double>& a, long b) {
  data += a.data * static_cast<double>(b);
}

template<>
inline void Z_NR<double>::submul(const Z_NR<double>& a, const Z_NR<double>& b) {
  data -= a.data * b.data;
}

template<>
inline void Z_NR<double>::submul_ui(const Z_NR<double>& a, unsigned long b) {
  data -= a.data * static_cast<double>(b);
}

template<>
inline void Z_NR<double>::abs(const Z_NR<double>& a) {
  data = fabs(a.data);
}

template<>
inline void Z_NR<double>::swap(Z_NR<double>& a) {
  std::swap(data, a.data);
}

/** random numbers */
template<>
inline void Z_NR<double>::randb(int bits) {
  mpz_t temp;
  mpz_init(temp);
  mpz_urandomb(temp, RandGen::getGMPState(), bits);
  data = mpz_get_d(temp);
  mpz_clear(temp);
}

template<>
inline void Z_NR<double>::randb_si(int bits) {
  randb (bits);
  data = data * RandGenInt::getBit();
}

template<>
inline void Z_NR<double>::randm(const Z_NR<double>& max) {
  mpz_t temp, lim;
  mpz_init(temp);
  mpz_init(lim);
  mpz_set_d(lim, max.data);
  mpz_urandomm(temp, RandGen::getGMPState(), lim);
  data = mpz_get_d(temp);
  mpz_clear(temp);
  mpz_clear(lim);
}

template<>
inline void Z_NR<double>::randm_si(const Z_NR<double>& max) {
  randm (max);
  data = data * RandGenInt::getBit();
}


template<>
inline void Z_NR<double>::nextprime(const Z_NR<double>& nbr) {
  mpz_t temp, temp2;
  mpz_init(temp);
  mpz_init(temp2);
  mpz_set_d(temp, nbr.data);
  mpz_nextprime(temp2, temp);
  data = mpz_get_d(temp2);
  mpz_clear(temp);
  mpz_clear(temp2);
}

/* operator Z_NR<double> */
template<>
inline ostream& operator<<(ostream& os, const Z_NR<double>& x) {
  return os << x.get_data();
}

template<>
inline istream& operator>>(istream& is, Z_NR<double>& x) {
  return is >> x.get_data();
}


#endif // #ifdef FPLLL_WITH_ZDOUBLE

FPLLL_END_NAMESPACE

#endif
