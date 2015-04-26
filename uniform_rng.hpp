#ifndef UNIFORM_RNG
#define UNIFORM_RNG
#include "distribution.hpp"
#include "basic_rng.hpp"
#include <cmath>
#include <cstdlib>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  class uniform_rng
    :public basic_rng<T_p,T_var>
  {
  private:
    T_var _min,_max;
  public:
    uniform_rng()
      :_min(0),_max(1)
    {}
    
    uniform_rng(const T_var& x1,
		const T_var& x2)
      :_min(x1),_max(x2)
    {}
    
    T_var do_get()const
    {
      return rand()/(T_var)RAND_MAX*(_max-_min)+_min;
    }

    basic_rng<T_p,T_var>* do_clone()const
    {
      return new uniform_rng(*this);
    }

    void do_var_range(T_var& x1,T_var& x2)const
    {
      x1=_min;
      x2=_max;
    }

    T_p do_eval(const T_var& x)const
    {
      if(x>=_min&&x<=_max)
	{
	  return 1;
	}
      else
	{
	  return 0;
	}
    }
  };
};

#endif
