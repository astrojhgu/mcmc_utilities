#ifndef DISCRETE_SAMPLE
#define DISCRETE_SAMPLE
#include "distribution.hpp"
#include "mcmc_exception.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

namespace mcmc_utilities
{

  template <typename T_p,typename T_var,typename T_urand>
  T_var discrete_sample(const probability_density_1d<T_p,T_var>& pd,const T_urand& urand)
  {
    std::vector<T_var> cp(pd.candidate_points());
    
    if(cp.empty())
      {
	throw no_candidate_points();
      }
    std::vector<T_p> prob(cp.size());
    prob[0]=std::exp(pd.eval_log(cp[0]));
    for(size_t i=1;i<cp.size();++i)
      {
	prob[i]=prob[i-1]+std::exp(pd.eval_log(cp[i]));
      }
    
    T_p p=urand()*prob.back();
    return cp[std::min<size_t>(std::distance(prob.begin(),std::upper_bound(prob.begin(),prob.end(),p)),cp.size()-1)];
    
  }
  
}

#endif
