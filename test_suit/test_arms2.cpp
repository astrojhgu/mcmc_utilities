#include <iostream>
#include <vector>
#include "core/gibbs_sampler.hpp"
using namespace std;
using namespace mcmc_utilities;

class std_norm
  :public probability_density_md<double,std::vector<double> >
{
  double do_eval_log(const std::vector<double>& x)const
  {
    //return exp(-(x*x/2.));
    double y=0;
    y=x[0]*x[0]+x[1]*x[1]+2*.9*x[1]*x[0];
    return (-y/(2*5*5.));
  }

  void do_var_range(double& x1,double& x2,const std::vector<double>& x,size_t ndim)const
  {
    x1=-100;
    x2=100;
  }

  
};

int main()
{
  std_norm sn;
  std::vector<double> factors(2,1);
  std::vector<double> init_var(2,0);
  srand(time(0));
  for(int i=0;i<10000;++i)
    {
      gibbs_sample(sn,init_var,0,u_random<double>);
      cout<<i<<" ";
      for(int j=0;j<init_var.size();++j)
	{
	  cout<<init_var[j]<<" ";
	}
      cout<<endl;
    }
}
