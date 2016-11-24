#ifndef GRAPH_ENSEMBLE_SAMPLE
#define GRAPH_ENSEMBLE_SAMPLE
#include "ensemble_sample.hpp"

namespace mcmc_utilities
{
  template <typename T,typename T_graph,typename T_ensemble>
  T_ensemble ensemble_sample(T_graph& g,
			      const T_ensemble& ensemble,
			      base_urand<T>& rng)
  {
    typedef typename element_type_trait<T_ensemble>::element_type T_var;
    auto p=g.get_params();
    size_t nparams=get_size(p);
    auto ensemble1=clone(ensemble);
    while(get_size(ensemble1)<2*nparams)
      {
	g.sample(rng);
	auto p=g.get_params();
	if(get_size(ensemble1)>0)
	  {
	    auto plast=get_element(ensemble1,get_size(ensemble1)-1);
	    bool no_jump=true;
	    for(size_t i=0;i<nparams;++i)
	      {
		if(get_element(p,i)!=get_element(plast,i))
		  {
		    no_jump=false;
		    break;
		  }
	      }
	    if(no_jump)
	      {
		continue;
	      }
	  }
	push_back(ensemble1,p);
      }
    ensemble1=ensemble_sample([&g](const T_var& x)
			     {
			       T result=g.eval_logprob(x);
			       return result;
			     },ensemble1,rng);
    return ensemble1;
  }
}


#endif
