/************************
 * Z=long specialization
 ************************/

#ifndef FPLLL_NR_Z_L_H
#define FPLLL_NR_Z_L_H

FPLLL_BEGIN_NAMESPACE

/* specialization long */
#ifdef FPLLL_WITH_ZLONG


template<>
inline Z_NR<long>::Z_NR() {}

template<>
inline Z_NR<long>::Z_NR(const Z_NR<long>& z) : data(z.data) {}

template<>
inline Z_NR<long>::~Z_NR() {}

/** get data */
template<>
inline double Z_NR<long>::get_d() const {
  return static_cast<double>(data);
}

#ifdef FPLLL_WITH_LONG_DOUBLE
template<>
inline long double Z_NR<long>::get_ld() const {
  return static_cast<long double>(data);
}
#endif

template<>
inline long Z_NR<long>::get_si() const {
  return data;
}

template<>
inline void Z_NR<long>::get_mpz(mpz_t r) const {
  mpz_set_si (r, data);
}

inline long compute_long_exponent(long data) {
  unsigned long y = static_cast<unsigned long>(abs(data));
  long e;
  for (e = 0; y; e++, y >>= 1) {}
  return e;
}

template<>
inline long Z_NR<long>::exponent() const {
  int int_expo;
  double f_norm = frexp(static_cast<double>(data), &int_expo);
  if (data > MAX_LONG_FAST && fabs(f_norm) == 0.5)
    return compute_long_exponent(data);
  else
    return static_cast<long>(int_expo);
}

/** set data */
template<>
inline void Z_NR<long>::set_str(const char* s) {
  data = atol(s);
}

/** comparison */
template<>
inline int Z_NR<long>::cmp(const Z_NR<long>& m) const {
  if (data > m.data) return 1;
  if (data == m.data) return 0;
  return -1;
}

template<>
inline int Z_NR<long>::sgn() const {
  if (data > 0) return 1;
  if (data == 0) return 0;
  return -1;
}

/** operator */
template<>
inline void Z_NR<long>::operator=(const Z_NR<long>& a) {
  data = a.data;
}

template<>
inline void Z_NR<long>::operator=(long a) {
  data = static_cast<long>(a);
}

template<>
inline void Z_NR<long>::operator=(const mpz_t& a) {
  data = mpz_get_si(a);
}

template<>
inline bool Z_NR<long>::operator<(const Z_NR<long>& a) const {
  return data < a.data;
}

template<>
inline bool Z_NR<long>::operator<(long a) const {
  return data < a;
}

template<>
inline bool Z_NR<long>::operator>(const Z_NR<long>& a) const {
  return data > a.data;
}

template<>
inline bool Z_NR<long>::operator>(long a) const {
  return data > a;
}

template<>
inline bool Z_NR<long>::operator<=(const Z_NR<long>& a) const {
  return data <= a.data;
}

template<>
inline bool Z_NR<long>::operator<=(long a) const {
  return data <= a;
}

template<>
inline bool Z_NR<long>::operator>=(const Z_NR<long>& a) const {
  return data >= a.data;
}

template<>
inline bool Z_NR<long>::operator>=(long a) const {
  return data >= a;
}

template<>
inline bool Z_NR<long>::operator==(const Z_NR<long>& a) const {
  return data == a.data;
}

template<>
inline bool Z_NR<long>::operator==(long a) const {
  return data == a;
}

template<>
inline bool Z_NR<long>::operator!=(const Z_NR<long>& a) const {
  return data != a.data;
}

template<>
inline bool Z_NR<long>::operator!=(long a) const {
  return data != a;
}

/** arithmetic */
template<>
inline void Z_NR<long>::add(const Z_NR<long>& a, const Z_NR<long>& b) {
  data = a.data + b.data;
}

template<>
inline void Z_NR<long>::add_ui(const Z_NR<long>& a, unsigned int b) {
  data = a.data + b;
}

template<>
inline void Z_NR<long>::sub(const Z_NR<long>& a, const Z_NR<long>& b) {
  data = a.data - b.data;
}

template<>
inline void Z_NR<long>::sub_ui(const Z_NR<long>& a, unsigned int b) {
  data = a.data - b;
}

template<>
inline void Z_NR<long>::neg(const Z_NR<long>& a) {
  data = -a.data;
}

template<>
inline void Z_NR<long>::mul(const Z_NR<long>& a, const Z_NR<long>& b) {
  data = a.data * b.data;
}

template<>
inline void Z_NR<long>::mul_si(const Z_NR<long>& a, long b) {
  data = a.data * b;
}

template<>
inline void Z_NR<long>::mul_ui(const Z_NR<long>& a, unsigned long b) {
  data = a.data * b;
}

template<>
inline void Z_NR<long>::mul_2si(const Z_NR<long>& a, long b) {
  //NOTE: if b >= size_in_bits(a), the result is undefined
  if (b >= 0)
    data = a.data << b;
  else
    data = a.data >> -b;
}

template<>
inline void Z_NR<long>::div_2si(const Z_NR<long>& a, long b) {
  //NOTE: if b >= size_in_bits(a), the result is undefined
  if (b >= 0)
    data = a.data >> b;
  else
    data = a.data << -b;
}

template<>
inline void Z_NR<long>::addmul(const Z_NR<long>& a, const Z_NR<long>& b) {
  data += a.data * b.data;
}

template<>
inline void Z_NR<long>::addmul_ui(const Z_NR<long>& a, unsigned long b) {
  data += a.data * b;
}

template<>
inline void Z_NR<long>::addmul_si(const Z_NR<long>& a, long b) {
  data += a.data * b;
}

template<>
inline void Z_NR<long>::submul(const Z_NR<long>& a, const Z_NR<long>& b) {
  data -= a.data * b.data;
}

template<>
inline void Z_NR<long>::submul_ui(const Z_NR<long>& a, unsigned long b) {
  data -= a.data * b;
}

template<>
inline void Z_NR<long>::abs(const Z_NR<long>& a) {
  data = labs(a.data);
}

template<>
inline void Z_NR<long>::swap(Z_NR<long>& a) {
  std::swap(data, a.data);
}

/** random numbers */
template<>
inline void Z_NR<long>::randb(int bits) {
  mpz_t temp;
  mpz_init(temp);
  mpz_urandomb(temp, RandGen::get_gmp_state(), bits);
  data = mpz_get_si(temp);
  mpz_clear(temp);
}

template<>
inline void Z_NR<long>::randb_si(int bits) {
  randb (bits);
  data = data * RandGenInt::get_bit();
}

template<>
inline void Z_NR<long>::randm(const Z_NR<long>& max) {
  mpz_t temp, lim;
  mpz_init(temp);
  mpz_init(lim);
  mpz_set_si(lim, max.data);
  mpz_urandomm(temp, RandGen::get_gmp_state(), lim);
  data = mpz_get_si(temp);
  mpz_clear(temp);
  mpz_clear(lim);
}

template<>
inline void Z_NR<long>::randm_si(const Z_NR<long>& max) {
  randm (max);
  data = data * RandGenInt::get_bit();
}


template<>
inline void Z_NR<long>::nextprime(const Z_NR<long>& nbr) {
  mpz_t temp, temp2;
  mpz_init(temp);
  mpz_init(temp2);
  mpz_set_ui(temp, nbr.data);
  mpz_nextprime(temp2, temp);
  data = mpz_get_ui(temp2);
  mpz_clear(temp);
  mpz_clear(temp2);
}


/* operators Z_NR<long> */
template<>
inline ostream& operator<<(ostream& os, const Z_NR<long>& x) {
  return os << x.get_data();
}

template<>
inline istream& operator>>(istream& is, Z_NR<long>& x) {
  return is >> x.get_data();
}


#endif // #ifdef FPLLL_WITH_ZLONG

FPLLL_END_NAMESPACE

#endif
