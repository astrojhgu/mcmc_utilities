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
    double k=param[0];
    double b=param[1];
    double s=param[2];
    for(int i=0;i<x_vec.size();++i)
      {
	std::vector<double> x(2);
	//log_p+=logdbivnorm(x,mu1,mu2,sigma11,sigma22,sigma12);
	log_p+=logdnorm(y_vec[i],k*x_vec[i]+b,s);
      }
    
    
    return log_p;
  }

  void do_var_range(double& x1,double& x2,const std::vector<double>& x0,size_t ndim)const
  {
    switch(ndim)
      {
      case 0:
	x1=-1;
	x2=1;
	break;
      case 1:
	x1=-1;
	x2=1;
	break;
      case 2:
	x1=.001;
	x2=10;
	break;
      default:
	assert(0);
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
  double xmin=1e99;
  double xmax=-1e99;
  ofstream ofs("fit_result.qdp");
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
      xmin=std::min(x,xmin);
      xmax=std::max(x,xmax);
      ofs<<x<<" "<<y<<endl;
    }
  //uniform_rng<double,double> ur;
  lm_model sn;
  std::vector<double> init_var(3,0);
  init_var[0]=1;
  init_var[1]=0;
  init_var[2]=.2;
  //cout<<sn.eval(init_var)<<endl;
  
  int s=(time(0));
  
  //cout<<"seed:"<<s<<endl;
  srand(s);
  u_random<double> rng;

  for(int i=0;i<10;++i)
    {
      gibbs_sample(sn,init_var,1,rng);
      cout<<i<<endl;
    }

  double kmean=0;
  
  for(int i=0;i<1000;++i)
    {
      gibbs_sample(sn,init_var,1,rng);
      
      kmean+=init_var[0];
      cout<<i<<" ";
      //cout<<init_var[1]<<" "<<init_var[2]<<endl;
      for(int j=0;j<init_var.size();++j)
	{
	  cout<<init_var[j]<<" ";
	}
      cout<<endl;
      cout.flush();
    }
  kmean/=1000;
  ofs<<"no no no"<<endl;
  for(int i=0;i<100;++i)
    {
      double x=(xmax-xmin)/100.*i+xmin;
      double y=kmean*x;
      ofs<<x<<" "<<y<<endl;
    }
  
  
}
