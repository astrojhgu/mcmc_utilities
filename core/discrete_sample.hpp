#ifndef DISCRETE_SAMPLE
#define DISCRETE_SAMPLE
#include "distribution.hpp"
#include "error_handler.hpp"
#include <vector>
#include <cmath>
#include "slicer.hpp"
#include <algorithm>

namespace mcmc_utilities
{

  template <typename T,typename TD,typename T_urand>
  T discrete_sample(const TD& pd,const std::pair<T,T>& xrange, const std::vector<T>& cp, T& xprev,size_t niter,T_urand& urand)
  {
    if(cp.empty())
      {
	//throw no_candidate_point();
	slice_sampler<T,TD> ss(pd,xrange,2,10,niter);
	xprev=ss.sample_double(xprev,urand);
	return xprev;
      }
    else
      {
	std::vector<T> prob(cp.size());
	prob[0]=std::exp(pd(cp[0]));
	for(size_t i=1;i<cp.size();++i)
	  {
	    prob[i]=prob[i-1]+std::exp(pd(cp[i]));
	  }
	
	T p=urand()*prob.back();
	return cp[std::min<size_t>(std::distance(prob.begin(),std::upper_bound(prob.begin(),prob.end(),p)),cp.size()-1)];
      }
  }
  
}

#endif
