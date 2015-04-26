#include "rejection_sampler_1d.hpp"
#include "uniform_rng.hpp"
#include <iostream>
#include "gibbs_sampler.hpp"
using namespace std;
using namespace mcmc_utilities;

class std_norm
  :public probability_density_md<double,double>
{
  double do_eval(const double& x)const
  {
    return exp(-(x*x/2.));
  }

  probability_density_md* do_clone()const
  {
    return new std_norm(*this);
  }

  void do_var_range(double& x1,double& x2)const
  {
    x1=-10;
    x2=10;
  }

  
};

int main()
{
  uniform_rng<double,double> ur;
  std_norm sn;
  
  for(int i=0;i<10000;++i)
    {
      cout<<rejection_sample(sn,uniform_rng<double,double>(-11,11),uniform_rng<double,double>(),40.)<<" "<<1<<endl;
    }
}
