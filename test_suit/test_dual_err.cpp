#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/distributions.hpp>
#include <uvsamplers/arms/arms_sampler.hpp>
#include <random/uniform.h>
#include <vector>
#include <set>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cassert>
#include <algorithm>
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
  int current_index;
  int size()const
  {
    return ndata;
  }
  variable(int n)
    :ndata(n*2+3),data(new double[ndata]),a(data[0]),b(data[1]),scat(data[2]),x(data+3),y(data+3+n),current_index(-1)
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

  bool relies_on(const double& v)const
  {
    if(&(data[current_index])==&v)
      {
	return true;
      }
    return false;
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

inline void set_element(variable& x,size_t i,const double& v)
{
  x[i]=v;
  x.current_index=i;
}


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
	//if(p.relies_on(p.x[i]))
	  {
	    log_p1+=logdnorm(x_vec[i],p.x[i],xe_vec[i]);
	  }
	//if(p.relies_on(p.y[i])||p.relies_on(p.b)||p.relies_on(p.a)||p.relies_on(p.x[i])||p.relies_on(p.scat))
	  {
	    log_p2+=logdnorm(p.y[i],p.b+p.a*(p.x[i]-2.3),p.scat);
	  }
	//if(p.relies_on(p.y[i]))
	  {
	    log_p3+=logdnorm(y_vec[i],p.y[i],ye_vec[i]);
	  }
      }
    log_p=log_p1+log_p2+log_p3;
    return log_p;
  }

  std::pair<double,double> do_var_range(const variable& x,size_t ndim)const
  {
    //for(int i=0;i<x_vec.size();++i)
    if(ndim>=3&&ndim<x_vec.size()+3)
      {
	return make_pair(-10,10);
      }
    else if(ndim>=x_vec.size()+3&&ndim<x_vec.size()*2+3)
      {
	return make_pair(-20,20);
      }
    else if(ndim==0)
      {
	return make_pair(0,20);
      }
    else if(ndim==1)
      {
	return make_pair(-20,20);
      }
    else
      {
	return make_pair(0.0001,100);
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
  arms_sampler<double,double> as;
  gibbs_sample(cd,p,as);
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,p,as);
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
