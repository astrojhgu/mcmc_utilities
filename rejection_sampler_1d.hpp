#ifndef REJ_SMP_MD
#define REJ_SMP_MD
#include "distribution.hpp"
#include "basic_rng.hpp"
#include "mcmc_exception.hpp"
#include <iostream>
namespace mcmc_utilities
{
  
  template <typename T_p,typename T_var>
  T_var rejection_sample(const probability_density_md<T_p,T_var>& pdensity,
		       const basic_rng<T_p,T_var>& enveloping_rng,
		       const basic_rng<T_p,T_var>& std_uniform_rng,
		       const T_p& factor)
  {
    T_var x1,x2;
    enveloping_rng.var_range(x1,x2);
    T_var y1,y2;
    pdensity.var_range(y1,y2);
    while(1)
      {
	T_var x=enveloping_rng.get();
	T_p pdensity_value=pdensity.eval(x);
	T_p enveloping_value=factor*enveloping_rng.eval(x);
	if(enveloping_value<pdensity_value)
	  {
	    std::cout<<"x="<<x<<" "<<enveloping_value<<" "<<pdensity_value<<std::endl;
	    throw not_enveloping();
	  }
	if(std_uniform_rng.get()*enveloping_value<pdensity_value)
	  {
	    return x;
	  }
      }
  };
  
}

#endif
