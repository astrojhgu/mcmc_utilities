#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/distributions.hpp>
#include <uvsamplers/arms/arms_sampler.hpp>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cassert>

using namespace std;
using namespace mcmc_utilities;
const double pi=4*atan(1);
struct variable
{
  typedef double value_type;
  double data[4];
  
  double& A;
  double& B;
  double& mu;
  double& sigma;

  int size()const
  {
    return 4;
  }
  
  variable()
    :A(data[0]),B(data[1]),mu(data[2]),sigma(data[3])
  {
    for(int i=0;i<3;++i)
      {
	data[i]=0;
      }
  }

  double& operator[](size_t i)
  {
    return data[i];
  }

  const double& operator[](size_t i)const
  {
    return data[i];
  }
};


class eff_distribution
  :public probability_density_md<double,variable>
{
private:
  std::vector<double> E;
  std::vector<int> nrec;
  std::vector<int> ninj;
public:
  eff_distribution()
  {
    ifstream ifs("eff.txt");
    for(;;)
      {
	double e1;
	int nrec1,ninj1;
	ifs>>e1>>nrec1>>ninj1;
	if(!ifs.good())
	  {
	    break;
	  }
	E.push_back(e1);
	nrec.push_back(nrec1);
	ninj.push_back(ninj1);
      }
  }

  double do_eval_log(const variable& x)const
  {
    double log_p=0;

    for(unsigned int i=0;i<E.size();++i)
      {
	double eff=x.A+(x.B-x.A)*phi((E[i]-x.mu)/x.sigma);
	//double eff=x.A+(x.B-x.A)*(atan((E[i]-x.mu)/x.sigma)+(pi/2))/pi;
	//dbin<double,double> d(eff,ninj[i]);
	double log_p1=logdbin(nrec[i],eff,ninj[i]);
	log_p+=log_p1;
      }
    return log_p;
  }

  std::pair<double,double> do_var_range(const variable& x,size_t ndim)const
  {
    if(ndim==0||ndim==1)
      {
	return make_pair(.001,1-1e-3);
      }
    else
      {
	return make_pair(.001,100);
      }
  }

  
  
};

int main()
{
  eff_distribution cd;
  variable x;
  
  x.A=.5;
  x.B=.75;
  x.mu=15;
  x.sigma=12;
  arms_sampler<double,double> as;
  //cout<<cd.eval_log(x)<<endl;
  for(int n=0;n<30000;++n)
    {
      gibbs_sample(cd,x,as);
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
