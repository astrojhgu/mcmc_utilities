#ifndef URAND_HPP
#define URAND_HPP
#include <cstdlib>
#include "base_urand.hpp"

namespace mcmc_utilities
{
  template <typename T>
  class urand
    :public base_urand<T>
  {
  private:
    T do_rand()const override
    {
      return rand()/(T)RAND_MAX;
    }
  };
}

#endif
