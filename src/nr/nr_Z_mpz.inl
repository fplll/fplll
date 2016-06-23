/****************************
 *  Z=mpz_t specialization
 ****************************/

#ifndef FPLLL_NR_Z_MPZ_H
#define FPLLL_NR_Z_MPZ_H

FPLLL_BEGIN_NAMESPACE


/* specialization */
template<>
inline Z_NR<mpz_t>::Z_NR() {
  mpz_init(data);
}

template<>
inline Z_NR<mpz_t>::Z_NR(const Z_NR<mpz_t>& z) {
  mpz_init_set(data, z.data);
}

template<>
inline Z_NR<mpz_t>::~Z_NR() {
  mpz_clear(data);
}

/** get data */
template<>
inline double Z_NR<mpz_t>::get_d() const {
  return mpz_get_d(data);
}

#ifdef FPLLL_WITH_LONG_DOUBLE
template<>
inline long double Z_NR<mpz_t>::get_ld() const {
  return LDConvHelper::mpz_get_ld(data);
}
#endif

template<>
inline long Z_NR<mpz_t>::get_si() const {
  return mpz_get_si(data);
}

template<>
inline void Z_NR<mpz_t>::get_mpz(mpz_t r) const {
  mpz_set (r, data);
}

template<>
inline long Z_NR<mpz_t>::exponent() const {
  long expo;
  mpz_get_d_2exp(&expo, data);
  return expo;
}

/** set data */
template<>
inline void Z_NR<mpz_t>::set_str(const char* s) {
  mpz_set_str(data, s, 10);
}

/** comparison */
template<>
inline int Z_NR<mpz_t>::cmp(const Z_NR<mpz_t>& m) const {
  int c = mpz_cmp(data, m.data);
  if (c > 0) return 1;
  if (c == 0) return 0;
  return -1;
}

template<>
inline int Z_NR<mpz_t>::sgn() const {
  return mpz_sgn(data);
}

/** operator */
template<>
inline void Z_NR<mpz_t>::operator=(const Z_NR<mpz_t>& a) {
  mpz_set(data, a.data);
}

template<>
inline void Z_NR<mpz_t>::operator=(const mpz_t& a) {
  mpz_set(data, a);
}

template<>
inline void Z_NR<mpz_t>::operator=(long a) {
  mpz_set_si(data, a);
}

template<>
inline bool Z_NR<mpz_t>::operator<(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) < 0;
}

template<>
inline bool Z_NR<mpz_t>::operator<(long a) const {
  return mpz_cmp_si(data, a) < 0;
}

template<>
inline bool Z_NR<mpz_t>::operator>(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) > 0;
}

template<>
inline bool Z_NR<mpz_t>::operator>(long a) const {
  return mpz_cmp_si(data, a) > 0;
}

template<>
inline bool Z_NR<mpz_t>::operator<=(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) <= 0;
}

template<>
inline bool Z_NR<mpz_t>::operator<=(long a) const {
  return mpz_cmp_si(data, a) <= 0;
}

template<>
inline bool Z_NR<mpz_t>::operator>=(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) >= 0;
}

template<>
inline bool Z_NR<mpz_t>::operator>=(long a) const {
  return mpz_cmp_si(data, a) >= 0;
}

template<>
inline bool Z_NR<mpz_t>::operator==(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) == 0;
}

template<>
inline bool Z_NR<mpz_t>::operator==(long a) const {
  return mpz_cmp_si(data, a) == 0;
}

template<>
inline bool Z_NR<mpz_t>::operator!=(const Z_NR<mpz_t>& a) const {
  return mpz_cmp(data, a.data) != 0;
}

template<>
inline bool Z_NR<mpz_t>::operator!=(long a) const {
  return mpz_cmp_si(data, a) != 0;
}

/** arithmetic */
template<>
inline void Z_NR<mpz_t>::add(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_add(data, a.data, b.data);
}

template<>
inline void Z_NR<mpz_t>::add_ui(const Z_NR<mpz_t>& a, unsigned int b) {
  mpz_add_ui(data, a.data, b);
}

template<>
inline void Z_NR<mpz_t>::sub(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_sub(data, a.data, b.data);
}

template<>
inline void Z_NR<mpz_t>::sub_ui(const Z_NR<mpz_t>& a, unsigned int b) {
  mpz_sub_ui(data, a.data, b);
}

template<>
inline void Z_NR<mpz_t>::neg(const Z_NR<mpz_t>& a) {
  mpz_neg(data, a.data);
}

template<>
inline void Z_NR<mpz_t>::mul(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_mul(data, a.data, b.data);
}

template<>
inline void Z_NR<mpz_t>::mul_si(const Z_NR<mpz_t>& a, long b) {
  mpz_mul_si(data, a.data, b);
}

template<>
inline void Z_NR<mpz_t>::mul_ui(const Z_NR<mpz_t>& a, unsigned long b) {
  mpz_mul_ui(data, a.data, b);
}

template<>
inline void Z_NR<mpz_t>::mul_2si(const Z_NR<mpz_t>& a, long b) {
  if (b >= 0)
    mpz_mul_2exp(data, a.data, b);
  else
    mpz_div_2exp(data, a.data, -b);
}

template<>
inline void Z_NR<mpz_t>::div_2si(const Z_NR<mpz_t>& a, long b) {
  if (b >= 0)
    mpz_div_2exp(data, a.data, b);
  else
    mpz_mul_2exp(data, a.data, -b);
}

template<>
inline void Z_NR<mpz_t>::addmul(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_addmul(data, a.data, b.data);
}

template<>
inline void Z_NR<mpz_t>::addmul_ui(const Z_NR<mpz_t>& a, unsigned long b) {
  mpz_addmul_ui(data, a.data, b);
}

template<>
inline void Z_NR<mpz_t>::addmul_si(const Z_NR<mpz_t>& a, long b) {
  if (b >= 0)
    mpz_addmul_ui(data, a.data, static_cast<unsigned long>(b));
  else
    mpz_submul_ui(data, a.data, static_cast<unsigned long>(-b));
}

template<>
inline void Z_NR<mpz_t>::submul(const Z_NR<mpz_t>& a, const Z_NR<mpz_t>& b) {
  mpz_submul(data, a.data, b.data);
}

template<>
inline void Z_NR<mpz_t>::submul_ui(const Z_NR<mpz_t>& a,unsigned long b) {
  mpz_submul_ui(data, a.data, b);
}

template<>
inline void Z_NR<mpz_t>::abs(const Z_NR<mpz_t>& x) {
  mpz_abs(data,x.data);
}

template<>
inline void Z_NR<mpz_t>::swap(Z_NR<mpz_t>& a) {
  mpz_swap(data, a.data);
}

/** random numbers */
template<>
inline void Z_NR<mpz_t>::randb(int bits) {
  mpz_urandomb(data, RandGen::getGMPState(), bits);
  if (bits > 32){ 
    unsigned long int tmp = mpz_get_ui(data);  
    gmp_randseed_ui(RandGen::gmpState, tmp*tmp);
    }
}

template<> 
inline void Z_NR<mpz_t>::randb_si(int bits)  {
  randb(bits);
  if (RandGenInt::getBit() < 0)
    mpz_neg(data, data);
}

template<>
inline void Z_NR<mpz_t>::randm(const Z_NR<mpz_t>& max) {
  mpz_urandomm(data, RandGen::getGMPState(), max.data);
}

template<>
inline void Z_NR<mpz_t>::randm_si(const Z_NR<mpz_t>& max) {
  randm(max);
  if (RandGenInt::getBit() < 0)
    mpz_neg(data, data);
}

template<>
inline void Z_NR<mpz_t>::nextprime(const Z_NR<mpz_t>& nbr) {
  mpz_nextprime(data, nbr.data);
}

/* FPLLL_V3_COMPAT */
#ifdef FPLLL_V3_COMPAT

template<>
inline void Z_NR<mpz_t>::print() const {
  mpz_out_str(stdout, 10, data);
}

template<>
inline void Z_NR<mpz_t>::printerr() const {
  mpz_out_str(stderr, 10, data);
}

template<>
inline void Z_NR<mpz_t>::read() {
  mpz_inp_str(data, stdin, 0);
}

template<>
inline double Z_NR<mpz_t>::get_d_2exp(long* expo) const {
  return mpz_get_d_2exp(expo,data);
}

template<>
inline void Z_NR<mpz_t>::set(/*const*/ mpz_t& d) {
  mpz_set(data,d);
}

template<>
inline void Z_NR<mpz_t>::set(unsigned long d) {
  mpz_set_ui(data,d);
}

#endif // #ifdef FPLLL_V3_COMPAT


/* operators Z_NR<mpz_t> */
template<>
inline ostream& operator<<(ostream& os, const Z_NR<mpz_t>& x) {
  int size = mpz_sizeinbase(x.getData(), 10) + 2;
  char* s = new char[size];
  mpz_get_str(s, 10, x.getData());
  os << s;
  delete [] s;
  return os;
}

template<>
inline istream& operator>>(istream& is, Z_NR<mpz_t>& x) {
  static string s;
  char c;

  s.clear();
  if (!(is >> c)) return is;
  if (c == '-') {
    s += c;
    if (!(is >> c)) return is;
  }
  if (!(c >= '0' && c <= '9')) {
    is.setstate(ios::failbit);
    return is;
  }

  do {
    s += c;
  } while (is.get(c) && c >= '0' && c <= '9');

  if (is.fail()) return is;
  if (is.eof())
    is.clear();
  else
    is.putback(c); 

  if (mpz_set_str(x.getData(), s.c_str(), 10)) {
    is.setstate(ios::failbit);
  }
  return is;
}


FPLLL_END_NAMESPACE

#endif
