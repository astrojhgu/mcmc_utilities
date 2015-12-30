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
  public:
    urand()
    {
      srand(time(0));
    }
  private:
    T do_rand()override
    {
      return rand()/(T)RAND_MAX;
    }
  };
}

#endif
