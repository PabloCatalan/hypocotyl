#ifndef _FUNCTIONS
#define _FUNCTIONS

#include <vector> //std::vector
#include <string> //std::string
#include <fstream> //std::ofstream
#include <map> //std::map
#include <cmath>//std::pow
#include <sstream>//std::stringstream
#include <random>
#include <chrono>

//BINOMIAL COEFFICIENT
/* double binomial_coefficient(int n, int k) { return 1/((n+1)*std::beta(n-k+1,k+1)); } */
template <class T = unsigned long>
T binomial_coefficient(unsigned long n, unsigned long k) {
    unsigned long i;
    T b;
    if (0 == k || n == k) {
        return 1;
    }
    if (k > n) {
        return 0;
    }
    if (k > (n - k)) {
        k = n - k;
    }
    if (1 == k) {
        return n;
    }
    b = 1;
    for (i = 1; i <= k; ++i) {
        b *= (n - (k - i));
        if (b < 0) return -1; /* Overflow */
        b /= i;
    }
    return b;
}
//INT_POW
int int_pow(int x, int p);
//INT_MIN
int int_min(int a, int b);
//EQUAL_VECTORS
bool equal_vectors(std::vector<int> vec1, std::vector<int> vec2);
//D_EQUAL
bool d_equal(double e1, double e2);
//D_LESS
bool d_less(double e1, double e2);
//D_MIN
double d_min(double e1, double e2);
//DECIMAL TO BINARY
inline std::string dectobin(int decimal,const int size){//converts decimal integers into binary std::strings
  std::string binary(size,'0');
  for (int i=size-1; i>-1; --i)
    if (decimal/std::pow(2,i)>=1){
      binary[size-i-1]='1';
      decimal=decimal-pow(2,i);
    }
  return binary;
}
//BINARY TO DECIMAL
inline int bintodec(const std::string& binary){//converts binary strings into decimal integers
  int decimal=0;
  int binary_size=binary.size();
  for(int i=0; i<binary_size; ++i)
    decimal += std::pow(2,binary_size-i-1)*(binary[i]-'0');
  return decimal;
}
//VECTOR TO STRING
inline std::string vec_to_str(const std::vector<int>& vec){
  std::string str(vec.size(),'0');
  for (int i=0; i<str.size(); ++i){
    std::stringstream ss;
    ss << vec[i];
    ss >> str[i];
  }
  return str;
}
//STRING TO VECTOR
inline std::vector<int> str_to_vec(const std::string& str){
  std::vector<int> vec(str.size(),0);
  for (int i=0; i<vec.size(); ++i)
    vec[i]=str[i]-'0';
  return vec;
}
//HAMMING
int vec_hamming(const std::vector<int>& v1, const std::vector<int>& v2);
int str_hamming(const std::string& v1, const std::string& v2);
//CHECK_FILE
void check_file(std::ifstream& file, const std::string& name);
//REVERSE
std::string reverse(const std::string& s1);
//IS IT DAY?
double is_day(double t, double Daylength);
//ELF3 PRODUCTION
double elf3p(double t, double Daylength, double pE1, double pE2, double L);
//GROWTH
std::vector<double> growth(const std::vector<double>& y, double t, int Temp, double Daylength, const std::string& mut, const std::vector<double>& params);
//RUNGE-KUTTA 4
std::vector<double> rk4(const std::vector<double>& u0, double t, double dt, int Temp, double Daylength, const std::string& mut, const std::vector<double>& params);
//MODEL RESULTS
std::pair<std::map<std::string, double>, std::map<std::string, std::vector<double> > > model_results(const std::vector<int>& Temp, const std::vector<double>& Daylength, const std::string& mut, const std::vector<double>& parameters);
//READ DATA
std::map<std::string, std::map<std::string, std::vector<double> > > read_data();
//READ EXPRESSION DATA
std::map<std::string, std::vector<double> > read_expression_data();
//ERROR
double error(const std::vector<double>& p, const std::map<std::string, std::map<std::string, std::vector<double> > >& data, const std::map<std::string, std::vector<double> >& expr, const std::vector<std::string>& mutants);
//CONSTRAINTS
void constraints(std::vector<double>& p2, std::default_random_engine& generator, std::uniform_real_distribution<double>& RNG, std::normal_distribution<double>& RNGnormal);
//SIMULATED ANNEALING
std::vector<double> simulated_annealing(int Iter, const std::string& suffix, std::vector<double>& parameters, const std::map<std::string, std::map<std::string, std::vector<double> > >& data, const std::map<std::string, std::vector<double> >& expr, const std::vector<std::string>& mutants, std::default_random_engine& generator, std::uniform_real_distribution<double>& RNG, std::normal_distribution<double>& RNGnormal, double& minerr);

#endif
