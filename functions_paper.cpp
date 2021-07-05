#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include "functions_paper.h"

const double pi = 3.14159265358979323846;

//INT_POW
int int_pow(int x, int p){
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = int_pow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}
//INT_MIN
int int_min(int a,int b){
  if (a<b)
    return a;
  else 
    return b;
}
//D_EQUAL
bool d_equal(double e1, double e2){//== operator for double
  return (std::abs(e1-e2)<1e-8);
}
//D_LESS
bool d_less(double e1, double e2){//< operator for double
  return (e1-e2<-1e-8);
}
//D_MIN
double d_min(double e1, double e2){
  if (d_less(e1,e2))
    return e1;
  else
    return e2;
}
//HAMMING
int vec_hamming(const std::vector<int>& v1, const std::vector<int>& v2){
  int ham=0;
  for (int i=0; i<v1.size(); ++i)
    if (v1[i]!=v2[i])
      ++ham;
  return ham;
}
int str_hamming(const std::string& v1, const std::string& v2){
  int ham=0;
  for (int i=0; i<v1.size(); ++i)
    if (v1[i]!=v2[i])
      ++ham;
  return ham;
}
//CHECK_FILE
void check_file(std::ifstream& file, const std::string& name){
  if (!file){
    std::cerr << "File " << name << " doesn't exist!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//IS IT DAY?
double is_day(double t, double Daylength){
  if (t<Daylength)
    return 1.0;
  else
    return 0.0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//ELF3 PRODUCTION
double elf3p(double t, double Daylength, double pE1, double pE2, double L){
  double k0=5;
  double t2=t-Daylength;
  double t3=t-24.0;
  if (d_equal(Daylength,0.0))
    return pE1+pE2;
  else if (d_equal(Daylength,24.0))
    return pE1-pE2;
  else
    return pE1-pE2*(-1.0+2.0/(1.0+std::exp(-k0*t))-2.0/(1.0+std::exp(-k0*t2))+2.0/(1.0+std::exp(-k0*t3)));
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//GROWTH
std::vector<double> growth(const std::vector<double>& y, double t, int Temp, double Daylength, const std::string& mut, const std::vector<double>& params){
  //PARAMETERS
  double pB28=params[0];
  double kr22=params[1];
  double kr28=params[2];
  double pE122=params[3];
  double pE128=params[4];
  double pE222=params[5];
  double pE228=params[6];
  double dE=params[7];
  double pPE22=params[8];
  double pPE28=params[9];
  double dP=params[10];
  double kPC=params[11];
  double dPB=params[12];
  double pCL28=params[13];
  double pCD=params[14];
  double dC=params[15];
  double pG=params[16];
  double kG=params[17];
  double pGP=params[18];
  double pGE=params[19];
  double pGB=params[20];
  double pGH=params[21];
  double pHC=params[22];
  //MUTPARS
  double mutBox=params[23];
  double mutEox=params[24];
  double mutPox=params[25];
  double mutPko1=params[26];
  double mutPko2=params[27];
  double mutCox=params[28];
  double mutCko1=params[29];
  double mutCko2=params[30];

  //VARIABLES
  double B=y[0];//PHYB
  double E=y[1];//ELF3
  double P=y[2];//PIF
  double C=y[3];//COP1
  double G=y[4];//Hypocotyl
  
  //ADAPT PARAMETERS
  double t1=std::fmod(t,24);
  double L=is_day(t1,Daylength);
  double pB=10.0;
  double kr=0.232;// datos de Casal
  double pE1=pE122;//adimensional
  double pE2=pE222;
  double pP=1.0;//adimensional
  double pPE=pPE22;
  double pCL=1.0;//adimensional
  double mB=1.0;//maximum PHYB value
  if (Temp==28){
    pB=pB28;
    kr=0.411;// datos de Casal
    pE1=pE128;
    pE2=pE228;
    pPE=pPE28;
    pCL=pCL28;
  }
  if (mut.find("PHYBox")!=std::string::npos)
    mB*=mutBox;
  if (mut.find("ELF3ox")!=std::string::npos)
    pE1*=mutEox;
  if (mut.find("PIF4ox")!=std::string::npos)
    pP*=mutPox;
  if (mut.find("pif4")!=std::string::npos)
    pP*=mutPko1;
  if (mut.find("pifq")!=std::string::npos)
    pP*=mutPko2;
  if (mut.find("COP1")!=std::string::npos){
    pCL*=mutCox;
    pCD*=mutCox;
  }
  if (mut.find("cop1-4")!=std::string::npos){
    pCL*=mutCko1;
    pCD*=mutCko1;
  }
  if (mut.find("cop1-6")!=std::string::npos){
    pCL*=mutCko2;
    pCD*=mutCko2;
  }
  if (mut.find("hy5")!=std::string::npos)
    pGH=0;   
  //Equations
  double dBdt=pB*L*(mB-B)-kr*B;
  double dEdt=elf3p(t1,Daylength,pE1,pE2,L)-dE*E;
  double dPdt=pP/(1+pPE*E)-dP*P/(1+kPC*C)-dPB*P*B;
  double dCdt=pCL*L+pCD*(1-L)-dC*C;
  double dGdt=pG+kG*pGP*P/(1+pGP*P+pGE*E+pGB*B+pGH/(1+pHC*C));
  std::vector<double> dydt={dBdt, dEdt, dPdt, dCdt, dGdt};
  if (mut.find("phyB")!=std::string::npos ||
      mut.find("phyb")!=std::string::npos)
    dydt[0]=0.0;        
  if (mut.find("elf3-8")!=std::string::npos)
    dydt[1]=0.0;        
  if (mut.find("cop1-null")!=std::string::npos)
    dydt[3]=0.0;     
  
  return dydt;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//RUNGE-KUTTA 4
std::vector<double> rk4(const std::vector<double>& u0, double t, double dt, int Temp, double Daylength, const std::string& mut, const std::vector<double>& params){
  //FIRST DERIVATIVE
  auto f0=growth(u0,t,Temp,Daylength,mut,params);
  //SECOND DERIVATIVE
  double t1=t+dt/2.0;
  std::vector<double> u1=u0;
  for (int i=0; i<u1.size(); ++i)
    u1[i]=u0[i]+dt*f0[i]/2.0;
  auto f1=growth(u1,t1,Temp,Daylength,mut,params);
  //THIRD DERIVATIVE
  double t2=t+dt/2.0;
  std::vector<double> u2=u0;
  for (int i=0; i<u2.size(); ++i)
    u2[i]=u0[i]+dt*f1[i]/2.0;
  auto f2=growth(u2,t2,Temp,Daylength,mut,params);
  //FOURTH DERIVATIVE
  double t3=t+dt;
  std::vector<double> u3=u0;
  for (int i=0; i<u3.size(); ++i)
    u3[i]=u0[i]+dt*f2[i];
  auto f3=growth(u3,t3,Temp,Daylength,mut,params);
  //COMBINE
  std::vector<double> u=u0;
  for (int i=0; i<u.size(); ++i)
    u[i]=u0[i]+dt*(f0[i]+2.0*f1[i]+2.0*f2[i]+f3[i])/6.0;
  return u;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//MODEL RESULTS
std::pair<std::map<std::string, double>, std::map<std::string, std::vector<double> > > model_results(const std::vector<int>& Temp, const std::vector<double>& Daylength, const std::string& mut, const std::vector<double>& parameters){
  std::map<std::string, double> rmap;
  std::map<std::string, double> rmap2;
  std::map<std::string,std::vector<double> > expr;
  std::vector<double> ttc(24);//times to check expression
  std::iota(ttc.begin(), ttc.end(), 48.0);
  double epsilon=1e-4;
  std::vector<std::string> mtc={"Col 22", "Col 28", "phyB 22", "phyB 28", "cop1-6 22", "cop1-6 28"};
  std::vector<double> x(5);
  for (int j=0; j<2; ++j){
    for (int l=0; l<5; ++l){
      //RESTART INITIAL CONDITIONS
      x[0]=0.0;//phyB
      x[1]=0.0;//ELF3
      x[2]=0.0;//PIF
      x[3]=0.0;//COP1
      x[4]=0.0;//HL
      //TEMPERATURE
      int T=Temp[j];
      std::stringstream sT;
      sT << T;
      std::string exprkey=mut+" "+sT.str();//to compare with expression data
      //LENGTH OF DAY
      double D=Daylength[l];
      std::stringstream sL;
      sL << D;
      //SYSTEM VARIABLES
      std::vector<std::vector<double> > x_vec;
      std::vector<double> times;
      x_vec.push_back(x);
      times.push_back(0.0);
      //SIMULATION
      int total_time=120;
      int datapoints=10;
      int Ttot=total_time*datapoints;
      double dt=1.0/(double)datapoints;
      double t=0.0;
      while (t<total_time){//for (int t=0; t<Ttot; ++t){
	std::vector<double> x=rk4(x_vec.back(),t,dt,T,D,mut,parameters);
	if (x[0]<0 || x[1]<0 || x[2]<0 || x[3]<0){
	  std::cout << "Something went wrong in " << mut << "\t" << T << "\t" << D << std::endl;
	  return std::make_pair(rmap2,expr);//concentrations of proteins can't be negative
	}
	double newtime=t+dt;
	//THIS IS FOR CHECKING EXPRESSION
	if (D==8 &&//it is short day
	    std::find(mtc.begin(), mtc.end(), exprkey)!=mtc.end() &&//it is a mutant to check
	    std::find_if(ttc.begin(), ttc.end(), [newtime](double b) { return std::abs(newtime - b) < 1e-8; })!=ttc.end()){//it is a time to check
	  expr[exprkey].push_back(x[1]);//ELF3 expression
	}
	x_vec.push_back(x);
	times.push_back(newtime);
	t+=dt;
      }
      std::string key=sT.str()+"_"+sL.str();
      rmap[key]=x_vec.back()[4];
    }//for daylengths
  }//for temperatures
  return std::make_pair(rmap, expr); 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//READ DATA
std::map<std::string, std::map<std::string, std::vector<double> > > read_data(){
  std::map<std::string, std::map<std::string, std::vector<double> > > data;
  std::ifstream in("data/daylength_def.csv");
  if (!in){
    std::cout << "No daylength.csv file!" << std::endl;
    std::exit(1);
  }
  std::vector<std::string> mutants;
  std::vector<std::string> ldata;
  std::vector<std::string> tdata;
  //READ DATA
  std::string line;
  std::getline(in, line);//names of columns
  std::string index, mut, temp, day;
  double growth;
  while (std::getline(in, line)){
    auto pos1=line.find(',');
    auto pos2=line.find(',', pos1+1);
    auto pos3=line.find(',', pos2+1);
    auto pos4=line.find(',', pos3+1);
    index=line.substr(0,pos1);
    mut=line.substr(pos1+1,pos2-pos1-1);
    day=line.substr(pos2+1,pos3-pos2-1);
    temp=line.substr(pos3+1,pos4-pos3-1);
    auto g=line.substr(pos4+1);
    std::stringstream ss;
    ss << g;
    ss >> growth;
    std::string key=temp+"_"+day;
    if (!data[mut].count(key))
      data[mut][key]=std::vector<double>();
    data[mut][key].push_back(growth);
  }
  in.close();
  return data;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//READ EXPRESSION DATA
std::map<std::string, std::vector<double> > read_expression_data(){
  std::map<std::string, std::vector<double> > data;
  std::ifstream in("data/expression_def.csv");
  if (!in){
    std::cout << "No expression.csv file!" << std::endl;
    std::exit(1);
  }
  std::vector<std::string> mutants;
  //READ KEYS
  std::string line;
  std::getline(in,line);
  std::stringstream ss;
  ss << line;
  std::string w, T;
  while (ss >> w >> T){
    std::string key=w+" "+T;
    mutants.push_back(key);
    data[key]=std::vector<double>();
  }
  //READ DATA POINTS
  while (std::getline(in,line)){
    std::stringstream ss;
    ss << line;
    double g1;
    int i=0;
    while (ss >> g1){
      data[mutants[i]].push_back(g1);
      ++i;
    }
  }
  in.close();
  return data;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//ERROR
double error(const std::vector<double>& p, const std::map<std::string, std::map<std::string, std::vector<double> > >& data, const std::map<std::string, std::vector<double> >& dataexpr, const std::vector<std::string>& mutants){
  //COMPUTE RESULTS AND COMPARE WITH DATA
  std::vector<double> length={0,8,12,16,24};
  std::vector<int> temp={22,28};
  double D=0.0;
  for (auto mut=mutants.begin(); mut!=mutants.end(); ++mut){
    auto results=model_results(temp, length, *mut, p);
    auto hypo=results.first;
    auto expr=results.second;
    if (hypo.empty()){
      std::cout << "something's wrong!" << std::endl;
      return 1e8;
    }
    //HYPOCOTYL DIFFERENCES
    for (auto it=hypo.begin(); it!=hypo.end(); ++it){
      auto key=it->first;
      auto pred=it->second;
      if (pred<0.0)
	return 1e8;
      if ((*mut).find("cop1-6")!=std::string::npos && key.back()=='0') continue;
      auto m1=data.at(*mut);
      if (!m1.count(key)) continue;
      int len=m1.at(key).size();
      for (auto it2=m1.at(key).begin(); it2!=m1.at(key).end(); ++it2){
	D+=(pred-*it2)*(pred-*it2)/(double)len;
      }
    }//for conditions
    //EXPRESSION DIFFERENCES
    double weight=1;
    for (auto it=expr.begin(); it!=expr.end(); ++it){
      auto key=it->first;
      auto pred=it->second;
      int ps=pred.size();
      auto m1=dataexpr.at(key);
      for (int t=0; t<ps; ++t){
    	D+=weight*(pred[t]-m1[t])*(pred[t]-m1[t]);
      }
    }
  }//for mutants
  return D;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//CONSTRAINTS
void constraints(std::vector<double>& p2, std::default_random_engine& generator, std::uniform_real_distribution<double>& RNG, std::normal_distribution<double>& RNGnormal){
  int psize=p2.size();
  int pos=RNG(generator)*psize;
  while (pos==1 || pos==2)//can't mutate kr
    pos=RNG(generator)*psize;
  //Choose movement
  double q=RNGnormal(generator);
  double try1=p2[pos]+q;
  while (try1<0.0){
    q=RNGnormal(generator);;
    try1=try1=p2[pos]+q;
  }
  p2[pos]=try1;
  //CONSTRAINTS
  // if (1.0<p2[5])//pE2 has to be less or equal pE1
  //   p2[5]=1.0;
  // if (p2[4]<p2[6])//pE2 has to be less or equal pE1
  //   p2[6]=p2[4];
  if (p2[23]<1.0)//mutBox has to be larger than 1
    p2[23]=1.01;
  if (p2[24]<1.0)//mutEox has to be larger than 1
    p2[24]=1.1;
  if (p2[25]<1.0)//mutPox has to be larger than 1
    p2[25]=1.1;
  if (p2[26]>1.0)//mutPko1 has to be smaller than 1
    p2[26]=0.9;
  if (p2[27]>1.0)//mutPko2 has to be smaller than 1
    p2[27]=0.9;
  if (p2[28]<1.0)//mutCox has to be larger than 1
    p2[28]=1.1;
  if (p2[29]>1.0)//mutCko1 has to be smaller than 1
    p2[29]=0.9;
  if(p2[30]>1.0)//mutCko2 has to be smaller than 1
    p2[30]=0.9;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//SIMULATED ANNEALING
std::vector<double> simulated_annealing(int Iter, const std::string& suffix, std::vector<double>& parameters, const std::map<std::string, std::map<std::string, std::vector<double> > >& data, const std::map<std::string, std::vector<double> >& expr, const std::vector<std::string>& mutants, std::default_random_engine& generator, std::uniform_real_distribution<double>& RNG, std::normal_distribution<double>& RNGnormal, double& minerror){
  double D=error(parameters, data, expr, mutants);//compute initial energy
  for (int i=0; i<Iter; ++i){//loop
    //Constraints on parameters
    std::vector<double> p2=parameters;
    constraints(p2, generator, RNG, RNGnormal);    
    //Compute new energy
    double D2=error(p2, data, expr, mutants);
    double ratio=D/D2;
    // double Ti=1.0/std::log(T0+i);
    double Ti=0.8/std::sqrt(i+1);
    std::cout << i << "\t" << D << "\t" << D2 << "\t" << ratio << "\t" << ratio*Ti << std::endl;
    if (ratio>=1.0){
      parameters=p2;
      D=D2;
    }
    else{
      double a1=RNG(generator);
      if (a1<ratio*Ti){
	parameters=p2;
	D=D2;
      }
    }
    //WRITE PARAMETERS EVERY TIME THEY REACH A MINIMU
    if (D<minerror){
      std::ofstream out("results/parameters_minerr_"+suffix+".txt");
      for (auto it=parameters.begin(); it!=parameters.end(); ++it)
	out << *it << std::endl;
      out.close();
    }
  }//for iterations
  return parameters;
}
