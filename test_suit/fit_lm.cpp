#include <iostream>
#include <vector>
#include <fstream>
#include <fio.h>
#include <math/distributions.hpp>
#include "core/gibbs_sampler.hpp"
using namespace std;
using namespace mcmc_utilities;
using namespace std;
std::vector<double> x_vec;
std::vector<double> y_vec;

class lm_model
  :public probability_density_md<double,std::vector<double> >
{
  double do_eval_log(const std::vector<double>& param)const
  {
    double log_p=0;
    double alpha=param[0];
    //double c=param[1];
    double c=10;
    double a=param[2];
    double b=param[3];
    double s=param[4];
    //log_p+=logdnorm(a,1.3,.5);
    //log_p+=logdnorm(b,0.,.5);
    //log_p+=logdnorm(s,.2,.5);
    //log_p+=logdnorm(c,1e1,1.);
    //log_p+=logdnorm(alpha,2.,1.);

    for(int i=0;i<x_vec.size();++i)
      {
	double x=x_vec[i];
	double y=y_vec[i];
	if(x<c)
	  {
	    assert(0);
	  }
	log_p+=logdpar(x,c,alpha);
	log_p+=logdlnorm(y,log(x)*a+b,s);
      }
    
    return log_p;
  }

  void do_var_range(double& x1,double& x2,const std::vector<double>& x0,size_t ndim)const
  {
    switch(ndim)
      {
      case 0:
	x1=1.;
	x2=5;
	break;
      case 1:
	x1=1;
	x2=10;
	break;
      case 2:
	x1=.5;
	x2=3.;
	break;
      case 3:
	x1=-1;
	x2=1;
	break;
      case 4:
	x1=.001;
	x2=1.;
	break;
      default:
	cerr<<ndim<<endl;
	assert(0);
	break;
      }
  }
};


int main(int argc,char* argv[])
{
  if(argc != 2)
    {
      return -1;
    }
  
  ifstream ifs(argv[1]);
  
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
  lm_model sn;
  std::vector<double> init_var(5,0);
  init_var[0]=4;
  init_var[1]=9.9;
  init_var[2]=1;
  init_var[3]=0;
  init_var[4]=.1;
  //cout<<sn.eval(init_var)<<endl;
  
  int s=(time(0));
  
  //cout<<"seed:"<<s<<endl;
  srand(s);
  for(int i=0;i<100000;++i)
    {

      cout<<i<<" ";
      //cout<<init_var[1]<<" "<<init_var[2]<<endl;
      for(int j=0;j<init_var.size();++j)
	{
	  cout<<init_var[j]<<" ";
	}
      cout<<endl;
      gibbs_sample(sn,init_var,1,u_random<double>);

    }
}
