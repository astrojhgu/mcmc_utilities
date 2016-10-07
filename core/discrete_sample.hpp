#ifndef DISCRETE_SAMPLE
#define DISCRETE_SAMPLE
#include "distribution.hpp"
#include "error_handler.hpp"
#include <vector>
#include <cmath>
#include "mcmc_traits.hpp"
#include "slicer.hpp"
#include <algorithm>

namespace mcmc_utilities
{

  template <typename T,typename TD,typename T_vector,typename T_urand>
  T discrete_sample(const TD& pd,const std::pair<T,T>& xrange, const T_vector& cp, T& xprev,size_t niter,T_urand& urand)
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
	T_vector prob(get_size(cp));
	set_element(prob,0,std::exp(pd(get_element(cp,0))));
	for(size_t i=1;i<get_size(cp);++i)
	  {
	    set_element(prob,i,get_element(prob,i-1)+std::exp(pd(get_element(cp,i))));
	  }
	
	T p=urand()*last_element(prob);
	return get_element(cp,std::min<size_t>(std::distance(std::begin(prob),std::upper_bound(begin(prob),std::end(prob),p)),get_size(cp)-1));
      }
  }
  
}

#endif
