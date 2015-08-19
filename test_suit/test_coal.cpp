#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/distributions.hpp>
#include <uvsamplers/arms/arms_sampler.hpp>
#include <vector>
#include <fstream>
#include <cassert>
using namespace std;
using namespace mcmc_utilities;

class coal_distribution
  :public probability_density_md<double,std::vector<double> >
{
private:
  std::vector<double> year;
  std::vector<double> n_accident;
public:
  coal_distribution()
  {
    ifstream ifs("coal.dat");
    for(;;)
      {
	double y,n;
	ifs>>y>>n;
	if(!ifs.good())
	  {
	    break;
	  }
	year.push_back(y);
	n_accident.push_back(n);
      }
  }

  double do_eval_log(const std::vector<double>& x)const
  {
    
    
    double lambda1=x[0];
    double lambda2=x[1];
    double k=x[2];
    double logp=0;
    double lambda=0;
    for(int i=0;i<year.size();++i)
      {
	if(year[i]<k)
	  {
	    lambda=lambda1;
	  }
	else
	  {
	    lambda=lambda2;
	  }
	int n=n_accident[i];
	double logp1=logdpoisson(n,lambda);
	logp+=logp1;
      }

    return logp;
  }

  std::pair<double,double> do_var_range(const std::vector<double>& x,size_t ndim)const
  {
    if(ndim==0||ndim==1)
      {
	return make_pair(0.0001,100.);
      }
    else
      {
	return make_pair(0,year.back()+1);
      }
  }
};

int main()
{
  coal_distribution cd;
  std::vector<double> x;
  x.push_back(3);
  x.push_back(4);
  x.push_back(50);
  //u_random<double> rng;
  arms_sampler<double,double> as;

  //cout<<cd.eval_log(x)<<endl;
  for(int n=0;n<10000;++n)
    {
      //gibbs_sample<double,std::vector<double> >(cd,x,1,as,10);
      gibbs_sample(cd,x,as);
      for(int i=0;i<x.size();++i)
	{
	  cout<<x[i]<<" ";
	}
      cout<<endl;
    }
}
