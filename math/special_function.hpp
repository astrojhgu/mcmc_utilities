#ifndef SPECIAL_FUNCTION_HPP
#define SPECIAL_FUNCTION_HPP
#include <cmath>
namespace mcmc_utilities
{
  
  template <typename T>
  T phi(T x)
  {
    static const T SQRT_2=std::sqrt((T)2.);
    return (1+erf(x/SQRT_2))/2;
  }
}


#endif
