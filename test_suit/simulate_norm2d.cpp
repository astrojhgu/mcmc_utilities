#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/distributions.hpp>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cassert>
using namespace std;

using namespace mcmc_utilities;
double mu1=0;
double mu2=0;
double sigma11=1;
double sigma22=1;
double sigma12=.3;

class norm2d
  :public probability_density_md<double,std::vector<double> >
{
public:
  norm2d()
  {
  }

  double do_eval_log(const std::vector<double>& x)const
  {
    double log_p=0;
    log_p+=logdbivnorm(x,mu1,mu2,sigma11,sigma22,sigma12);
    return log_p;
  }

  //void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  void do_var_range(double& x1,double& x2,const std::vector<double>& x,size_t ndim)const
  {
    x1=-1e2;
    x2=1e2;
  }
};

int main()
{
  norm2d cd;
  std::vector<double> var;
  var.push_back(1);
  var.push_back(1);
  u_random<double> rng;
  double xmin,xmax;

  for(int n=0;n<1000;++n)
    {
      gibbs_sample(cd,var,1,rng,100); 
    }
  for(long n=0;n<1000000l;++n)
    {
      gibbs_sample(cd,var,1,rng,100);
      if(n%100==0)
	{
	  for(unsigned int i=0;i<var.size();++i)
	    {
	      //cout<<std::log(var[i])<<" ";
	      cout<<var[i]<<" ";
	    }
	  cout<<endl;
	}
    }
}
