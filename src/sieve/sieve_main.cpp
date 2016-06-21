#include "sieve_main.h"

long dot_time;
long dot_num;
long count_bad;

/**
 * help function
 */
static void
main_usage (char *myself) {
  cout << "Usage: " << myself << " [options]\n"
       << "List of options:\n"
       << "  -a [2|3|4]\n"
       << "     2- or 3- or 4-sieve;\n"
       << "  -f filename\n"
       << "     Input filename\n"
       << "  -r nnn\n"
       << "     Generate a random instance of dimension nnn\n"
       << "  -t nnn\n"
       << "     Targeted norm^2=nnn\n"
       << "  -v\n"
       << "     Verbose mode\n";
  exit(0);
}


/**
 * run sieve
 */
template<class ZT>
int main_run_sieve (ZZ_mat<ZT> B, Z_NR<ZT> goal_norm, int alg)
{
  Gauss_sieve<ZT, FP_NR<double> > gsieve(B, alg);
  gsieve.set_goal_norm2 (goal_norm);
  if (gsieve.alg == 3)
    gsieve.run_3sieve();
  else if (gsieve.alg == 4)
    gsieve.run_4sieve();
  else
    gsieve.run_2sieve();

  return 0;
}


/**
 * main function
 */
int main (int argc, char** argv)
{
  char* input_file_name = NULL;
  char* goal_norm_s = NULL;
  bool flag_verbose = false, flag_file = false;
  int option, alg, dim = 10;

  dot_time = 0;
  dot_num = 0;
  count_bad = 0;
  alg = 2;

  /* parse */
  while((option = getopt (argc, argv, "a:f:r:t:v")) != -1) {
    switch (option) {
    case 'a':
      alg = atoi(optarg);
      if (alg != 2 && alg != 3 && alg != 4) {
        cout << " Error, only support 2-, 3- and 4-sieve" <<endl;
        exit(1);
      }
    case 'f':
      input_file_name = optarg;
      flag_file = true;
      break;
    case 'r':
      dim = atoi(optarg);
      flag_file = false;
      break;
    case 'v':
      flag_verbose = true;
      break;
    case 't':
      //ngoal_norm = atol(optarg);
      cout << optarg << endl;
      goal_norm_s = optarg;
      break;
    case 'h':
    case '?':
      main_usage(argv[0]);
      return -1;
    }
  }

  /* set lattice */
  ZZ_mat<mpz_t> B(dim,dim);
  if (flag_file) {
    ifstream input_file(input_file_name);
    if (input_file.is_open()) {
      input_file >> B;
      input_file.close();
    }
    else {
      cin >> B;
    }
    cout << "# [info] reading lattice of dimension "
         << B.getRows() << "x" << B.getCols() << endl;
  }
  else {
    cout << "# [info] generating random lattice of dimension "
         << dim << endl;
    //srand (time(NULL));
    srand (1);
    B.gen_intrel(10*dim);
  }

  /* set targeted norm */
  Z_NR<mpz_t> goal_norm, max;
  if (goal_norm_s != NULL) {
    goal_norm.set_str(goal_norm_s);
    //goal_norm.mul(goal_norm, goal_norm);
  }
  if (goal_norm < 0)
    goal_norm = 0;
  cout << "# [info] goal norm^2 is " << goal_norm << endl;
  
  /* preprocessing of basis */
  clock_t stime = clock();
  lllReduction(B, 0.75, 0.51, LM_WRAPPER);
  //bkzReduction(B, 20, BKZ_DEFAULT, FT_DEFAULT, 0);

  clock_t etime = clock();
  double secs = (etime - stime) / (double) CLOCKS_PER_SEC;
  cout << "# [info] LLL took time " << secs << " s" << endl;
  //cout << B << endl;

  /* decide integer type */
  stime = clock();
  max = B.getMax();

#if 1
  if (max < std::numeric_limits<int>::max()) {
    long goal_norm_l = abs(goal_norm.get_si());
    Z_NR<long> goal_norm_lt;
    goal_norm_lt = goal_norm_l;
    ZZ_mat<long> B2(B.getRows(),B.getCols());
    for (int i = 0; i < B.getRows(); i ++)
      for (int j = 0; j < B.getCols(); j ++)
        B2(i, j) = B(i, j).get_si();
    main_run_sieve<long>(B2, goal_norm_lt, alg);
  }
  else
#endif
    main_run_sieve<mpz_t>(B, goal_norm, alg);
  
  etime = clock();
  secs = (etime - stime) / (double) CLOCKS_PER_SEC;
  cout << "# [info] sieve took time " << secs << " s" << endl;

  /* dot product time */
  cout << "# [info] dot_time " << dot_time << endl;
  cout << "# [info] dot_num " << dot_num << endl;
  cout << "# [info] dot_time/dot_number " << (double) dot_time/dot_num << endl;

  return 1;
}
