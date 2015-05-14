#ifndef BASIC_RNG_HPP
#define BASIC_RNG_HPP
#include "distribution.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  class basic_rng
    :public probability_density_md<T_p,T_var>
  {
  public:
    T_var get()const
    {
      return do_get();
    }

  private:
    virtual T_var do_get()const=0;
  };
};
#endif

