#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/distributions.hpp>
#include <random/uniform.h>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cassert>

using namespace std;
using namespace mcmc_utilities;

struct variable
{
  const int ndata;
  double* data;
  double& a;
  double& b;
  double& scat;
  typedef double value_type;
  double *x,*y;
  int size()const
  {
    return ndata;
  }
  variable(int n)
    :ndata(n*2+3),data(new double[ndata]),a(data[0]),b(data[1]),scat(data[2]),x(data+3),y(data+3+n)
  {}

  double& operator[](size_t i)
  {
    assert(i<ndata);
    return data[i];
  }

  const double& operator[](size_t i)const
  {
    assert(i<ndata);
    return data[i];
  }
public:
  ~variable()
  {
    delete[] data;
  }
private:
  variable(const variable& rhs);
  variable operator=(const variable& rhs);
};


class dual_err_distribution
  :public probability_density_md<double,variable >
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

  double do_eval_log(const variable& p)const
  {
    double log_p=0;
    
    double intrscat=1/std::sqrt(p.scat);
    double log_prior_a=logdt(p.a,0.,1.,1.);
    double log_prior_b=logdnorm(p.b,0.,100.);
    double log_prior_scat=logdgamma(p.scat,1e-2,1e-2);

    log_p+=log_prior_a;
    log_p+=log_prior_b;
    log_p+=log_prior_scat;
    double log_p1=0,log_p2=0,log_p3=0;
    for(int i=0;i<x_vec.size();++i)
      {
	log_p1+=logdnorm(x_vec[i],p.x[i],xe_vec[i]);
	log_p2+=logdnorm(p.y[i],p.b+p.a*(p.x[i]-2.3),p.scat);
	log_p3+=logdnorm(y_vec[i],p.y[i],ye_vec[i]);
      }
    log_p=log_p1+log_p2+log_p3;
    return log_p;
  }

  void do_var_range(double& x1,double& x2,const variable& x,size_t ndim)const
  {
    //for(int i=0;i<x_vec.size();++i)
    if(ndim>=3&&ndim<x_vec.size()+3)
      {
	x1=-10;x2=10;
      }
    else if(ndim>=x_vec.size()+3&&ndim<x_vec.size()*2+3)
      {
	x1=-20;
	x2=20;
      }
    else if(ndim==0)
      {
	x1=0;
	x2=20;
      }
    else if(ndim==1)
      {
	x1=-20;
	x2=20;
      }
    else
      {
	x1=0;
	x2=100;
      }
  }  
};

int main()
{
  dual_err_distribution cd;
  variable p(cd.x_vec.size());
  for(int i=0;i<cd.x_vec.size();++i)
    {
      //x.push_back(cd.y_vec[i]);
      p.x[i]=cd.x_vec[i];
      p.y[i]=cd.y_vec[i];
    }
  p.a=3;
  p.b=8;
  p.scat=.2;
  //gibbs_sample(cd,x);
  ranlib::Uniform<double> uf;
  //cout<<cd.eval_log(x)<<endl;
  gibbs_sample(cd,p,1,[&uf](){return uf.random();},2);
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,p,1,[&uf](){return uf.random();},10);
      if(n>100)
	{
	  for(int i=0;i<3;++i)
	    {
	      cout<<p[i]<<" ";
	    }
	  cout<<endl;
	}
    }
}
