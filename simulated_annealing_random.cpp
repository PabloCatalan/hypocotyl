#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <map>
#include <chrono>
#include <random>
#include "functions_paper.h"

//MAIN
int main(int argc, char* argv[]){
  
  //RANDOM NUMBER GENERATOR
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  
  //WE START THE PROCESS argv[1] TIMES
  std::stringstream ssloop;
  ssloop << argv[1];
  int Loop;
  ssloop >> Loop;

  //LENGTH OF LOOP
  std::stringstream ssiter;
  ssiter << argv[2];
  int Iter;
  ssiter >> Iter;
  
  //SUFFIX
  std::string suffix;
  std::stringstream ssuf;
  ssuf << argv[3];
  ssuf >> suffix;
  
  //READ DATA
  std::map<std::string, std::map<std::string, std::vector<double> > > data=read_data();
  std::map<std::string, std::vector<double> > expr=read_expression_data();

  //PARAMETERS
  std::normal_distribution<double> RNGs (1.0,0.1);
  std::vector<double> parameters(31);
  for (int i=0; i<parameters.size(); ++i){
    double r1=RNGs(generator);
    while (r1<0)
      r1=RNGs(generator);
    parameters[i]=r1;
  }
  parameters[1]=0.232;
  parameters[2]=0.411;
  //CHECK IF OTHER PARAMETERS EXIST
  std::ifstream in("results/parameters_minerr_"+suffix+".txt");
  if (in){
    parameters.clear();
    double p1;
    while (in >> p1)
      parameters.push_back(p1);
  }
  in.close();    
  //MUTANTS
  std::vector<std::string> mutants;
  for (auto it=data.begin(); it!=data.end(); ++it){
    if (it->first=="det1-1") continue;
    //else if (it->first=="hy5") continue;
    else if (it->first=="elf4") continue;
    else if (it->first=="lux") continue;
    mutants.push_back(it->first);
  }
  //RANDOM NUMBER GENERATOR
  std::uniform_real_distribution<double> RNG (0.0,1.0);
  //MINIMAL ERROR
  double minerr=error(parameters, data, expr, mutants);//compute initial energy
  //LOOP
  for (int k=0; k<Loop; ++k){
    std::cout << k << std::endl;
    std::normal_distribution<double> RNGnormal(0.0,0.05);
    //SIMULATED ANNEALING
    std::vector<double> p2=simulated_annealing(Iter, suffix, parameters, data, expr, mutants, generator, RNG, RNGnormal, minerr);
    std::ofstream out("results/parameters_"+suffix+"_random_noELF3exp.txt");
    for (auto it=parameters.begin(); it!=parameters.end(); ++it)
      out << *it << std::endl;
    out.close();
    parameters=p2;
  }
}
