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
  //SUFFIX
  std::string suffix;
  std::stringstream ssuf;
  ssuf << argv[1];
  ssuf >> suffix;
  
  //READ DATA
  std::map<std::string, std::map<std::string, std::vector<double> > > data=read_data();
  std::map<std::string, std::vector<double> > expr=read_expression_data();

  //PARAMETERS
  std::string ifile="results/parameters_"+suffix+".txt";
  std::vector<double> parameters;
  {
    std::ifstream in(ifile);
    check_file(in, ifile);
    double p1;
    while (in >> p1)
      parameters.push_back(p1);
  }
  //MUTANTS
  std::vector<std::string> mutants;
  for (auto it=data.begin(); it!=data.end(); ++it){
    if (it->first=="det1-1") continue;
    //else if (it->first=="hy5") continue;
    else if (it->first=="elf4") continue;
    else if (it->first=="lux") continue;
    mutants.push_back(it->first);
  }
  double D=error(parameters, data, expr, mutants);
  std::ofstream out("results/error_"+suffix+".txt");
  out << D << std::endl;
  std::cout << "El error es " << D << std::endl;
  out.close();
}
