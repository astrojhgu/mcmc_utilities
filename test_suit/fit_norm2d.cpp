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
    double mu1=param[0];
    double mu2=param[1];
    double sigma11=param[2];
    double sigma22=param[3];
    double sigma12=param[4];

    for(int i=0;i<x_vec.size();++i)
      {
	std::vector<double> x(2);
	x[0]=x_vec[i];
	x[1]=y_vec[i];
	log_p+=logdbivnorm(x,mu1,mu2,sigma11,sigma22,sigma12);
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
      case 3:
	x1=.001;
	x2=10;
	break;
      case 4:
	x1=-10;
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
  init_var[0]=0;
  init_var[1]=0;
  init_var[2]=1;
  init_var[3]=1;
  init_var[4]=0;
  //cout<<sn.eval(init_var)<<endl;
  
  int s=(time(0));
  
  //cout<<"seed:"<<s<<endl;
  srand(s);
  u_random<double> rng;
  for(int i=0;i<100000;++i)
    {

      cout<<i<<" ";
      //cout<<init_var[1]<<" "<<init_var[2]<<endl;
      for(int j=0;j<init_var.size();++j)
	{
	  cout<<init_var[j]<<" ";
	}
      cout<<endl;
      gibbs_sample(sn,init_var,1,rng);

    }
}
