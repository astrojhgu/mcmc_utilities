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
double a=1.5;
double b=0.2;
double s=.1;
double c=10;
double alpha=3.5;

class LM
  :public probability_density_md<double,std::vector<double> >
{
public:
  LM()
  {
  }

  double do_eval_log(const std::vector<double>& x)const
  {
    double log_p=0;
    log_p+=logdpar(x[0],c,alpha);
    //log_p+=logdnorm(x[0],500.,50.);
    //log_p+=logdnorm(std::log(x[1]),log(x[0])*a+b,s);
    log_p+=logdlnorm(x[1],log(x[0])*a+b,s);
    return log_p;
  }

  //void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  void do_var_range(double& x1,double& x2,const std::vector<double>& x,size_t ndim)const
  {
    x1=1e1;
    x2=1e6;
  }
};

int main()
{
  LM cd;
  std::vector<double> var;
  var.push_back(4000);
  var.push_back(exp(log(var[0])*a+b));
  u_random<double> rng;
  size_t idx=1;
  double xmin,xmax;

  for(int n=0;n<1000;++n)
    {
      gibbs_sample(cd,var,1,rng,100); 
    }
  for(long n=0;n<300000l;++n)
    {
      for(int i=0;i<var.size();++i)
	{
	  gibbs_sample(cd,var,i,true,rng,10,true);
	}
      //if(n%100==0)
	{
	  for(unsigned int i=0;i<var.size();++i)
	    {
	      //cout<<std::log(var[i])<<" ";
	      cout<<var[i]<<" ";
	    }
	  cout<<endl;
	  //cout<<var[0]<<" "<<var[1]<<"\n";
	  xmax=max(xmax,var[0]);
	  cout.flush();
	}
    }
  return 0;
  cout<<"no no no"<<endl;

  for(double x=1E1;x<xmax;x*=1.1)
    {
      cout<<x<<" "<<exp(a*log(x)+b)<<endl;
    }
}
