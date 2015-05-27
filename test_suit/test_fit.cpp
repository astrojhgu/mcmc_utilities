#include <iostream>
#include <vector>
#include <fstream>
#include <fio.h>
#include "core/gibbs_sampler.hpp"
using namespace std;
using namespace mcmc_utilities;
using namespace std;
std::vector<double> x_vec;
std::vector<double> y_vec;

class variable
  :public std::vector<double>
{
public:
  double& k;
  double& b;
  double& sigma;

  variable()
    :vector(3),k(at(0)),b(at(1)),sigma(at(2))
  {}

private:
  variable(const variable&);
  variable& operator=(const variable&);
};

class prob
  :public probability_density_md<double,variable>
{
  double do_eval_log(const variable& param)const
  {
    double p=0;
    
    for(int i=0;i<x_vec.size();++i)
      {
	double y1=param.k*x_vec[i]+param.b;
	p+=(-(y_vec[i]-y1)*(y_vec[i]-y1)/(2*param.sigma*param.sigma)-std::log(param.sigma));
      }
    
    return p;
  }

  void do_var_range(double& x1,double& x2,const variable& x0,size_t ndim)const
  {
    x2=20;
    if(ndim==2)
      {
	x1=.001;
      }
    else
      {
	x1=-20;
      }
  }

  
};


int main()
{
  ifstream ifs("data.txt");
  for(;;)
    {
      double x,y;
      ifs>>x>>y;
      if(!ifs.good())
	{
	  break;
	}
      x_vec.push_back(x);
      y_vec.push_back(y);
    }
  //uniform_rng<double,double> ur;
  prob sn;
  variable init_var;
  init_var[0]=4;
  init_var[1]=3;
  init_var[2]=20;
  
  //cout<<sn.eval(init_var)<<endl;
  
  int s=(time(0));
  
  //cout<<"seed:"<<s<<endl;
  srand(s);
  //return 0;
  init_var[0]=4;
  init_var[1]=3;
  init_var[2]=.5;
  u_random<double> rng;
  for(int i=0;i<100000;++i)
    {
      gibbs_sample(sn,init_var,1,rng);

      cout<<i<<" ";
      //cout<<init_var[1]<<" "<<init_var[2]<<endl;
      for(int j=0;j<init_var.size();++j)
	{
	  cout<<init_var[j]<<" ";
	}
      cout<<endl;
    }
}
