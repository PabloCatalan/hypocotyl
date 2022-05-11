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
    mutants.push_back(it->first);
  }
  //SENSITIVITY
  double dS=0.00005;
  double D0=error(parameters, data, expr, mutants);
  std::ofstream out("results/sensitivity_analysis_global_"+suffix+".txt");
  for (int p=0; p<parameters.size(); ++p){
    //PERTURBATION 1
    std::vector<double> p2=parameters;
    p2[p]+=dS*p2[p];
    double D1=error(p2, data, expr, mutants);
    //PERTURBATION 2
    std::vector<double> p3=parameters;
    p3[p]-=dS*p3[p];
    double D2=error(p3, data, expr, mutants);
    //DERIVATIVE
    double dY=(D1-D2)/D0/(2*dS*parameters[p]);
    std::cout << p << "\t" << D1 << "\t" << D2 << "\t" << D0 << "\t" << dY << std::endl;
    if (dY>1e5)
      out << p << "," << std::endl;
    else
      out << p << "," << dY << std::endl;
  }
  out.close();
}
