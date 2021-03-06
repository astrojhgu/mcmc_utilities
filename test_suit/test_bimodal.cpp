#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/distributions.hpp>
#include <core/ensemble_sample.hpp>
#include <core/ptsample.hpp>
#include <random>
#include <vector>
#include <fstream>
#include <cassert>
using namespace std;
using namespace mcmc_utilities;

template <typename T>
using std_vector=std::vector<T>;
int main()
{
  srand(time(0));
  std::random_device rd;
  std::mt19937 gen(rd());
  
  
  std::vector<std::vector<std::vector<double> > > ensemble_list;
  constexpr int nwalker=100;
  constexpr int ntemp=10;
  std::vector<double> beta_list;
  for(int i=0;i<ntemp;++i)
    {
      if(i==ntemp-1)
	{
	  beta_list.push_back(0);
	}
      else
	{
	  beta_list.push_back(pow(1/std::sqrt(2),(double)i));
	}
      
      std::vector<std::vector<double> > ensemble;
      for(int j=0;j<nwalker;++j)
	{
	  std::vector<double> x{urng<double>(gen)*20-10,urng<double>(gen)};
	  ensemble.push_back(x);
	}
      ensemble_list.push_back(ensemble);
    }
  
      
  
  //cout<<cd.eval_log(x)<<endl;
  std::vector<std::vector<double> > results(ntemp);
  
  for(int n=0;n<100000;++n)
    {
      if(n%100==0)
	{
	  std::cerr<<n<<std::endl;
	}
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
	  double mu1=-5;
	  double mu2=5;
	  double sigma1=.1;
	  double sigma2=1;
	  
	  if(x[1]<.5){return -(x[0]-mu1)*(x[0]-mu1)/(2*sigma1*sigma1)-std::log(sigma1);}
	  else{return -(x[0]-mu2)*(x[0]-mu2)/(2*sigma2*sigma2)-std::log(sigma2);};
	},ensemble_list,gen,beta_list,n%10==0,1);
      int j=0;
      do{j=urng<double>(gen)*nwalker;}while(j>=nwalker);
      //if(n%100==0)
	{
	  for(int i=0;i<ntemp;++i)
	    {
	      results[i].push_back(ensemble_list[i][0][0]);
	    }
	}

    }
  for(auto& i:results)
    {
      std::sort(i.begin(),i.end());
      for(int j=0;j<i.size();++j)
	{
	  std::cout<<i[j]<<" "<<j<<std::endl;
	}
      std::cout<<"no no no"<<std::endl;
    }

}
