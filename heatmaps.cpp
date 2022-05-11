#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <map>
#include "functions_paper.h"

//MAIN
int main(int argc, char* argv[]){

  //MUTANT
  std::vector<std::string> M={"COP1ox WT", "COP1ox phyb", "COP1ox elf3-8", 
			      "PHYBox WT", "PHYBox cop1-4", "PHYBox elf3-8",
			      "ELF3ox WT", "ELF3ox phyb", "ELF3ox cop1-4"};
  
  //SUFFIX
  std::string suffix;
  std::stringstream ssuf;
  ssuf << argv[1];
  ssuf >> suffix;  

  //INITIAL CONDITIONS
  std::vector<double> x(5);
  x[0]=0.0;//PHYB
  x[1]=0.0;//ELF3
  x[2]=0.0;//PIF4
  x[3]=0.0;//COP1
  x[4]=0.0;//Hypocotyl

  //PARAMETERS
  std::vector<double> parameters;
  std::ifstream in("results/parameters_"+suffix+".txt");
  check_file(in, "results/parameters_"+suffix+".txt");
  double p1;
  while (in >> p1)
    parameters.push_back(p1);

  //LOOP OVER ALL TEMPERATURES
  std::vector<int> temp(2);
  temp[0]=22;
  temp[1]=28;
  std::vector<double> length;
  double step1=0.25;//0.25 si queremos más resolución
  for (double d=0; d<24.5; d=d+step1)
    length.push_back(d);
  std::vector<double> prot;//{0.01, 0.1, 1.0, 10.0, 100.0, 1000.0};
  double step2=0.05;//0.05 si queremos más resolución
  for (double c=-3; c<3; c=c+step2)
    prot.push_back(std::pow(10,c));
  //LOOP OVER ALL MUTANTS
  for (auto it=M.begin(); it!=M.end(); ++it){
    std::string mut=*it;
    std::cout << mut << std::endl;
    std::ofstream out("results/heatmap_"+mut+"_"+suffix+".csv");
    out << "Prot,T,D,Growth" << std::endl;
    for (int c=0; c<prot.size(); ++c){
      for (int j=0; j<temp.size(); ++j){
	for (int l=0; l<length.size(); ++l){
	  auto p2=parameters;
	  //RESTART INITIAL CONDITIONS
	  x[0]=0.0;//phyB
	  x[1]=0.0;//ELF3
	  x[2]=0.0;//PIF
	  x[3]=0.0;//COP1
	  x[4]=0.0;//HL
	  //TEMPERATURE
	  int T=temp[j];
	  std::stringstream sT;
	  sT << T;
	  //LENGTH OF DAY
	  double D=length[l];
	  std::stringstream sL;
	  sL << D;
	  if (mut.find("COP1ox")!=std::string::npos)//COP1
	    p2[28]=prot[c];
	  else if (mut.find("PHYBox")!=std::string::npos)//phyB
	    p2[23]=prot[c];
	  else if (mut.find("ELF3ox")!=std::string::npos)//ELF3
	    p2[24]=prot[c];
	  else if (mut.find("PIF4ox")!=std::string::npos)//PIF4
	    p2[25]=prot[c];
	  //SYSTEM VARIABLES
	  std::vector<std::vector<double> > x_vec;
	  std::vector<double> times;
	  x_vec.push_back(x);
	  times.push_back(0.0);
	  //SIMULATION
	  int total_time=120;
	  int datapoints=100;
	  int Ttot=total_time*datapoints;
	  double dt=1.0/(double)datapoints;
	  for (int t=0; t<Ttot; ++t){
	    double time=t*dt;
	    std::vector<double> x=rk4(x_vec[t],time,dt,T,D,mut,p2);
	    x_vec.push_back(x);
	    times.push_back(time+dt);
	  }
	  out << prot[c] << "," << T << "," << D << "," << x_vec.back()[4] << std::endl;
	}//for daylengths
      }//for temperatures
    }//for prot
    out.close();
  }//for mutants  
}
