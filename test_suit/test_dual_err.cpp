#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/special_function.hpp>
#include <random/uniform.h>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cassert>

using namespace std;
using namespace mcmc_utilities;

class dual_err_distribution
  :public probability_density_md<double,std::vector<double> >
{
public:
  std::vector<double> x_vec;
  std::vector<double> xe_vec;
  std::vector<double> y_vec;
  std::vector<double> ye_vec;
public:
  dual_err_distribution()
  {
    ifstream ifs("dual_err.txt");
    for(;;)
      {
	double x,y,xe,ye;
	ifs>>x>>xe>>y>>ye;
	if(!ifs.good())
	  {
	    break;
	  }
	x_vec.push_back(x);
	xe_vec.push_back(xe);
	y_vec.push_back(y);
	ye_vec.push_back(ye);
      }
  }

  double do_eval_log(const std::vector<double>& x)const
  {
    double log_p=0;
    double a=x[0];
    double b=x[1];
    double scat=x[2];
    double intrscat=1/std::sqrt(scat);
    double log_prior_a=logdt(a,0.,1.,1.);
    double log_prior_b=logdnorm(b,0.,100.);
    double log_prior_scat=logdgamma(scat,1e-2,1e-2);

    log_p+=log_prior_a;
    log_p+=log_prior_b;
    log_p+=log_prior_scat;
    double log_p1=0,log_p2=0,log_p3=0;
    for(int i=0;i<x_vec.size();++i)
      {
	log_p1+=logdnorm(x_vec[i],x[i+3],xe_vec[i]);
	log_p2+=logdnorm(x[3+i+x_vec.size()],b+a*(x[i+3]-2.3),scat);
	log_p3+=logdnorm(y_vec[i],x[i+3+x_vec.size()],ye_vec[i]);
      }
    log_p=log_p1+log_p2+log_p3;
    return log_p;
  }

  void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  {
    x1.resize(x_vec.size()*2+3);
    x2.resize(x_vec.size()*2+3);
    for(int i=0;i<x_vec.size();++i)
      {
	x1[i+3]=-10;x2[i+3]=10;
	x1[i+3+x_vec.size()]=-20;
	x2[i+3+x_vec.size()]=20;
	
      }
    x1[0]=0;
    x2[0]=10;

    x1[1]=-10;
    x2[1]=10;

    x1[2]=0;
    x2[2]=100;
  }  
};

int main()
{
  dual_err_distribution cd;
  std::vector<double> x(cd.x_vec.size()*2+3);
  for(int i=0;i<cd.x_vec.size();++i)
    {
      //x.push_back(cd.y_vec[i]);
      x[i+3]=cd.x_vec[i];
      x[i+3+cd.x_vec.size()]=cd.y_vec[i];
    }
  x[0]=3;
  x[1]=8;
  x[2]=.2;
  //gibbs_sample(cd,x);
  ranlib::Uniform<double> uf;
  //cout<<cd.eval_log(x)<<endl;
  gibbs_sample(cd,x,1,[&uf](){return uf.random();},2);
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,x,1,[&uf](){return uf.random();},10);
      if(n>100)
	{
	  for(int i=0;i<3;++i)
	    {
	      cout<<x[i]<<" ";
	    }
	  cout<<endl;
	}
    }
}
