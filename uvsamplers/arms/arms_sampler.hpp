#ifndef ARMSSAMPLER_HPP
#define ARMSSAMPLER_HPP
#include <core/uvsampler.hpp>
#include <arms.hpp>
#include <cstdlib>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  class arms_sampler
    :public uvsampler<T_p,T_var>
  {
    class{
    public:
      T_p operator()()const
      {
	return rand()/(T_var)RAND_MAX;
      }
    }rnd;
    
    T_var do_sample(const probability_density_1d<T_p,T_var>& pd,
		    const T_var& xprev)const override
    {
      std::vector<T_var> xsamp(10);
      arms_simple(pd,xprev,xsamp,true,rnd);
      //return arms(pd,xprev,10,rnd);
      return xsamp.back();
    }
  };
};


#endif

