#ifndef CONTINUOUS_SAMPLE
#define CONTINUOUS_SAMPLE

#include "distribution.hpp"
#include "error_handler.hpp"
#include <vector>
#include <cmath>
#include "slicer.hpp"
#include "arms.hpp"
#include <algorithm>

namespace mcmc_utilities
{
  template <typename T,typename T_urand>
  T continuous_sample(const probability_density_1d<T>& pd,T& xprev,size_t niter,const T_urand& urand)
  {
    size_t xmchange_count=0;
    
    xprev=arms(pd,xprev,niter,urand,xmchange_count);
    
    if(xmchange_count==0)
      {
	slice_sampler<T> ss(pd,1,niter,10);
	xprev=ss.sample_step(xprev,urand);
      }
    
    return xprev;
  }
}


#endif
