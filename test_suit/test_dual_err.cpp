#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/special_function.hpp>
#include <distribution/dbin.hpp>
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
    double a=x[x_vec.size()];
    double b=x[x_vec.size()+1];
    double scat=x[x_vec.size()+2];
    double intrscat=1/std::sqrt(scat);
    double log_prior_a=logdt(a,0.,1.,1.);
    double log_prior_b=logdnorm(b,0.,1e-4);
    double log_prior_scat=logdgamma(scat,1e-2,1e-2);

    log_p+=log_prior_a;
    log_p+=log_prior_b;
    log_p+=log_prior_scat;

    for(int i=0;i<x_vec.size();++i)
      {
	log_p+=logdnorm(x_vec[i],x[i],std::pow(xe_vec[i],-2.));
	log_p+=logdnorm(y_vec[i],b+a*(x[i]-2.3),1/(ye_vec[i]*ye_vec[i]+scat*scat));
      }
    return log_p;
  }

  dual_err_distribution* do_clone()const
  {
    return new dual_err_distribution(*this);
  }

  void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  {
    x1.resize(x_vec.size()+3);
    x2.resize(x_vec.size()+3);
    for(int i=0;i<x_vec.size();++i)
      {
	x1[i]=-10;x2[i]=10;
      }
    x1[x_vec.size()]=0;
    x2[x_vec.size()]=10;

    x1[x_vec.size()+1]=-10;
    x2[x_vec.size()+1]=10;

    x1[x_vec.size()+2]=0;
    x2[x_vec.size()+2]=1;
  }  
};

int main()
{
  dual_err_distribution cd;
  std::vector<double> x(cd.x_vec);
  x.push_back(4);
  x.push_back(8.2);
  x.push_back(.3);

  //cout<<cd.eval_log(x)<<endl;
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,x);
      //if(n>100)
	{
	  for(unsigned int i=0;i<x.size();++i)
	    {
	      cout<<x[i]<<" ";
	    }
	  cout<<endl;
	}
    }
}
