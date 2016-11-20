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
int main()
{
  urand<double> rng;
  prng<double> prng;
  
  std::vector<std::vector<std::vector<double> > > ensemble_list;
  constexpr int nwalker=100;
  constexpr int ntemp=10;
  std::vector<double> beta_list;
  for(int i=0;i<ntemp;++i)
    {
      beta_list.push_back(pow(.5,(double)i));
      std::vector<std::vector<double> > ensemble;
      for(int j=0;j<nwalker;++j)
	{
	  std::vector<double> x{rng()*20-10,rng()};
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
      
      ensemble_list=ptsample([](const std::vector<double>& x){
	  if(x[0]<-15||x[0]>15||x[1]<0||x[1]>1)
	    {
	      return -numeric_limits<double>::infinity();
	    }
	  
	  if(x[1]<.5){return -(x[0]+5)*(x[0]+5)/(2*.1*.1);}
	  else{return -(x[0]-5)*(x[0]-5);};
	},ensemble_list,prng,beta_list,n%2==0);
      int j=0;
      do{j=rng()*nwalker;}while(j>=nwalker);
      //if(n%100==0)
	{
	  for(int i=0;i<1;++i)
	    {
	      cout<<ensemble_list[0][0][i]<<" ";
	    }
	  cout<<endl;
	}
    }
}
