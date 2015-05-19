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
double a=1.3;
double b=0;
double c=.2;
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
    log_p+=-2*log(x[0])+11;
    log_p+=logdnorm(std::log10(x[1]),log10(x[0])*a+b,c);
    return log_p;
  }

  //void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  void do_var_range(double& x1,double& x2,const std::vector<double>& x,size_t ndim)const
  {
    //x1[0]=1E1;x2[0]=1E5;
    //x1[1]=1E1;x2[1]=1E5;
    //x1=1e1;
    //x2=1e5;
    static double xmin[2],xmax[2];
    xmin[0]=1E1;
    xmin[1]=pow(10.,log10(xmin[0])*a+b-5*c);
    xmax[0]=1E7;
    xmax[1]=pow(10.,log10(xmax[0])*a+b+5*c);
    if(ndim==0)
      {
	x1=pow(10.,(log10(x[1])-b)/a-5*c);
	x2=pow(10.,(log10(x[1])-b)/a+5*c);
	if(x1>x2)
	  {
	    double y=x1;
	    x1=x2;
	    x2=y;
	  }
	x1=max(x1,xmin[1]);
	x2=min(x2,xmax[1]);
      }
    else
      {
	x1=pow(10.,log10(x[0])*a+b-5*c);
	x2=pow(10.,log10(x[0])*a+b+5*c);
	if(x1>x2)
	  {
	    double y=x1;
	    x1=x2;
	    x2=y;
	  }
	x1=max(x1,xmin[0]);
	x2=min(x2,xmax[0]);
	
      }
  }

  
  
};

int main()
{
  LM cd;
  std::vector<double> x;
  x.push_back(2e2);
  x.push_back(pow(10.,log10(x[0])*a+b));
  //cout<<std::log(std::numeric_limits<float>::max()/10)<<endl;
  //exit(0);
  //cout<<cd.eval_log(x)<<endl;
  double xmax=0;
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,x,1,u_random<double>,1000);
      //if(n>100)
	{
	  for(unsigned int i=0;i<x.size();++i)
	    {
	      //cout<<std::log10(x[i])<<" ";
	    }
	  cout<<x[0]<<" "<<x[1]<<"\n";
	  xmax=max(xmax,x[0]);
	  cout.flush();
	}
    }

  cout<<"no no no"<<endl;

  for(double x=1E1;x<xmax;x*=1.1)
    {
      cout<<x<<" "<<pow(10.,a*log10(x)+b)<<endl;
    }
}
