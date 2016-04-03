#ifndef CONTINUOUS_SAMPLE
#define CONTINUOUS_SAMPLE

//#include "distribution.hpp"
#include "error_handler.hpp"
#include <vector>
#include <cmath>
#include "slicer.hpp"
#include "arms.hpp"
#include <algorithm>

namespace mcmc_utilities
{
  template <typename T,typename TD, typename T_urand>
  T continuous_sample(const TD& pd,const std::pair<T,T>& xrange, const std::vector<T>& init_x, T& xprev,size_t niter,T_urand& urand)
  {
    size_t xmchange_count=0;
    try
      {
	xprev=arms(pd,xrange, init_x, xprev,niter,urand,xmchange_count);
      }
    catch(const ill_conditioned_distribution& e)
      {
	xmchange_count=0;
      }
    catch(const too_many_rejections& e)
      {
	xmchange_count=0;
      }
    if(xmchange_count==0)
      {
	slice_sampler<T,TD> ss(pd,xrange, 1,niter,10);
	xprev=ss.sample_step(xprev,urand);
      }
    
    return xprev;
  }
}


#endif
