#ifndef PRNG_HPP
#define PRNG_HPP
#include "prng_engine.hpp"

namespace mcmc_utilities
{
  template <typename T>
  class prng
    :public base_urand<T>
  {
  private:
    sitmo::prng_engine eng;
  public:
    prng()
    {}

    prng(int n)
      :eng(n)
    {}

  private:
    T do_rand()const override
    {
      return const_cast<sitmo::prng_engine&>(eng)()/static_cast<T>(sitmo::prng_engine::max());
    }
    
  };
}


#endif
