#include <iostream>
#include <vector>
#include "gibbs_sampler.hpp"
using namespace std;
using namespace mcmc_utilities;

class std_norm
  :public probability_density_md<double,std::vector<double> >
{
  double do_eval(const std::vector<double>& x)const
  {
    //return exp(-(x*x/2.));
    double y=0;
    y=x[0]*x[0]+x[1]*x[1]+2*.9*x[1]*x[0];
    return exp(-y/(2*5*5.));
  }

  probability_density_md* do_clone()const
  {
    return new std_norm(*this);
  }

  void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  {
    for(int i=0;i<x1.size();++i)
      {
	x1[i]=-100;
	x2[i]=100;
      }
  }

  
};

int main()
{
  std_norm sn;
  std::vector<double> factors(2,1);
  std::vector<double> init_var(2,0);
  std::vector<double> xmin,xmax;
  sn.var_range(xmin,xmax);
  srand(time(0));
  for(int i=0;i<10000;++i)
    {
      gibbs_sample(sn,init_var,0);
      cout<<i<<" ";
      for(int j=0;j<init_var.size();++j)
	{
	  cout<<init_var[j]<<" ";
	}
      cout<<endl;
    }
}
