#ifndef NORMAL_RNG_HPP
#define NORMAL_RNG_HPP

#include <random>
#include <cmath>
#include "../core/base_urand.hpp"

namespace mcmc_utilities
{
  template <typename T>
  T normal_rng(base_urand<T>& urand)
  {
    static std::normal_distribution<T> dist((T)0,(T)1);
    return dist(urand);
  }

  template <typename T>
  T normal_rng(T m,T s,base_urand<T>& urand)
  {
    return normal_rng(urand)*s+m;
  }
};



#endif
