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

class dugongs
  :public probability_density_md<double,std::vector<double> >
{
private:
  std::vector<double> x_vec;
  std::vector<double> y_vec;
public:
  dugongs()
  {
    ifstream ifs("dugongs.dat");
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
  }

  double do_eval_log(const std::vector<double>& p)const
  {
    double alpha=p[0];
    double beta=p[1];
    double U3=p[2];
    double tau=p[3];
    double sigma=1/sqrt(tau);
    //double U3=logit(gamma);
    double gamma=ilogit(U3);
    double log_p=0;

    log_p+=logdnorm(alpha,0.,pow(1e-6,-2.));
    log_p+=logdnorm(beta,0.,pow(1e-6,-2.));
    log_p+=logdgamma(tau,.001,.001);
    log_p+=logdnorm(U3,0.,pow(1e-4,-2.));
    for(int i=0;i<x_vec.size();++i)
      {
	double mu=alpha-beta*pow(gamma,x_vec[i]);
	log_p+=logdnorm(y_vec[i],mu,tau);
      }
    return log_p;
  }

  void do_var_range(double& x1,double& x2,const std::vector<double>& x,size_t ndim)const
  {
    switch(ndim)
      {
      case 0:
      case 1:
	x1=-5;
	x2=5;
	break;
      case 2:
	x1=.1;
	x2=4;
	break;
      case 3:
	x1=1e-10;
	x2=200;
	break;
      default:
	assert(0);
      }
  }  
};

int main()
{
  dugongs cd;
  std::vector<double> x;
  x.push_back(2.6);
  x.push_back(.9);
  x.push_back(.87);
  x.push_back(100);
  u_random<double> rng;
  //cout<<cd.eval_log(x)<<endl;
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,x,1,rng);
      if(n>1000)
	{
	  for(unsigned int i=0;i<x.size();++i)
	    {
	      cout<<x[i]<<" ";
	    }
	  cout<<endl;
	}
    }
}
