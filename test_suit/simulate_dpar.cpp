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
double c=10;
double alpha=2;
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
    
    log_p+=logdpar(x[0],c,alpha);
    return log_p;
  }

  //void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  void do_var_range(double& x1,double& x2,const std::vector<double>& x,size_t ndim)const
  {
    x1=10;
    x2=1e5;
  }

  
  
};

int main()
{
  LM cd;
  std::vector<double> x;
  x.push_back(2e2);
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
	      cout<<x[i]<<" ";
	    }
	  cout<<endl;
	  //cout<<x[0]<<" "<<x[1]<<"\n";
	}
    }

}
