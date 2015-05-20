#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/special_function.hpp>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cassert>
using namespace std;

using namespace mcmc_utilities;
double a=1;
double b=0;
double s=.4;
double c=10;
double alpha=4;

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
    //log_p+=g*log(x[0]);
    //assert(x[0]>0);
    log_p+=logdpar(x[0],c,alpha);
    //log_p+=logdnorm(std::log(x[1]),log(x[0])*a+b,s);
    log_p+=logdlnorm(x[1],log(x[0])*a+b,s);
    return log_p;
  }

  //void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  void do_var_range(double& x1,double& x2,const std::vector<double>& x,size_t ndim)const
  {
    x1=1e1;
    x2=1e6;
    return;
  }

  
  
};

int main()
{
  LM cd;
  std::vector<double> x;
  x.push_back(2e1);
  x.push_back(exp(log(x[0])*a+b));
  //cout<<std::log(std::numeric_limits<float>::max()/10)<<endl;
  //exit(0);
  //cout<<cd.eval_log(x)<<endl;
  double xmax=0;
  for(int n=0;n<100000;++n)
    {
      gibbs_sample(cd,x,1,u_random<double>,100);
      //if(n>100)
	{
	  for(unsigned int i=0;i<x.size();++i)
	    {
	      //cout<<std::log(x[i])<<" ";
	      cout<<x[i]<<" ";
	    }
	  cout<<endl;
	  //cout<<x[0]<<" "<<x[1]<<"\n";
	  xmax=max(xmax,x[0]);
	  cout.flush();
	}
    }

  cout<<"no no no"<<endl;

  for(double x=1E1;x<xmax;x*=1.1)
    {
      cout<<x<<" "<<exp(a*log(x)+b)<<endl;
    }
}
