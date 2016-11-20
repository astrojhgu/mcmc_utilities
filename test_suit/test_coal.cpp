#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/distributions.hpp>
#include <core/ensemble_sample.hpp>
#include <core/ptsample.hpp>
#include <core/urand.hpp>
#include <rng/prng.hpp>
#include <vector>
#include <fstream>
#include <cassert>
using namespace std;
using namespace mcmc_utilities;

template <typename T>
using std_vector=std::vector<T>;


class coal_distribution
  :public probability_density_md<double,std::vector<double>,std_vector>
{
private:
  std::vector<double> year;
  std::vector<double> n_accident;
public:
  coal_distribution()
    :year(),n_accident()
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

  double do_eval_log(const std::vector<double>& x,int n)const
  {   
    double lambda1=x[0];
    double lambda2=x[1];
    double k=x[2];

    if(lambda1<=0||lambda2<=0||k<0||k>100)
      {
	return -std::numeric_limits<double>::infinity();
      }
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
    //logp+=200;
    //std::cout<<logp<<std::endl;
    return logp;
  }

  std::pair<double,double> do_var_range(const std::vector<double>& x,size_t ndim)const
  {
    if(ndim==0||ndim==1)
      {
	return make_pair(1e-6,100.);
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
  urand<double> rng;
  prng<double> prng;
  
  std::vector<std::vector<std::vector<double> > > ensemble_list;
  constexpr int nwalker=100;
  constexpr int ntemp=10;
  std::vector<double> beta_list;
  for(int i=0;i<ntemp;++i)
    {
      beta_list.push_back(i/(ntemp-1));
      std::vector<std::vector<double> > ensemble;
      for(int j=0;j<nwalker;++j)
	{
	  std::vector<double> x{6+rng()-.5,1+rng()-.5,60+rng()*4-2};
	  ensemble.push_back(x);
	}
      ensemble_list.push_back(ensemble);
    }
  
      
  
  //cout<<cd.eval_log(x)<<endl;
  for(int n=0;n<10000;++n)
    {
      //gibbs_sample<double,std::vector<double> >(cd,x,1,as,10);
      //gibbs_sample(cd,x,rng);
      /*
      ensemble=ensemble_sample([&cd](const std::vector<double>& x){
	  double y=cd.eval_log(x);
	  //std::cerr<<x[0]<<" "<<y<<std::endl;
	  return y;},ensemble,prng,4);
      */
      
      ensemble_list=ptsample([&cd](const std::vector<double>& x){
	  double y=cd.eval_log(x);
	  //std::cerr<<x[0]<<" "<<y<<std::endl;
	  return y;},ensemble_list,prng,beta_list,n%10==0);
      int j=0;
      do{j=rng()*nwalker;}while(j>=nwalker);
      for(int i=0;i<3;++i)
	{
	  cout<<ensemble_list.back()[j][i]<<" ";
	}
      cout<<endl;
    }
}
