#ifndef DISCRETE_SAMPLE
#define DISCRETE_SAMPLE
#include "distribution.hpp"
#include "mcmc_exception.hpp"
#include <vector>
#include <cmath>
#include "slicer.hpp"
#include <algorithm>

namespace mcmc_utilities
{

  template <typename T,typename T_urand>
  T discrete_sample(const probability_density_1d<T>& pd,T& xprev,const T_urand& urand)
  {
    std::vector<T> cp(pd.candidate_points());
    
    if(cp.empty())
      {
	//throw no_candidate_point();
	slice_sampler<T> ss(pd,2,10);
	xprev=ss.sample(xprev,urand);
	return xprev;
      }
    else
      {
	std::vector<T> prob(cp.size());
	prob[0]=std::exp(pd.eval_log(cp[0]));
	for(size_t i=1;i<cp.size();++i)
	  {
	    prob[i]=prob[i-1]+std::exp(pd.eval_log(cp[i]));
	  }
	
	T p=urand()*prob.back();
	return cp[std::min<size_t>(std::distance(prob.begin(),std::upper_bound(prob.begin(),prob.end(),p)),cp.size()-1)];
      }
  }
  
}

#endif
