#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP
#include <cmath>
#include "../core/mcmc_exception.hpp"
namespace mcmc_utilities
{
  
  template <typename T>
  T phi(T x)
  {
    static const T SQRT_2=std::sqrt((T)2.);
    return (1+erf(x/SQRT_2))/2;
  }

  template <typename T>
  T lbeta(T x,T y)
  {
    return lgamma(x)+lgamma(y)-lgamma(x+y);
  }

  template <typename T_p,typename T_var>
  T_p log_factorial(const T_var& n)
  {
    //http://en.wikipedia.org/wiki/Factorial
    return lgamma(n+1);
  }

  template <typename T_p,typename T_var>
  T_p log_CN(const T_var& m,const T_var& n)
  {
    T_p result= log_factorial<T_p,T_var>(m)-log_factorial<T_p,T_var>(n)-log_factorial<T_p,T_var>(m-n);
    return result;
  }
  
  template <typename T>
  T logit(const T& x)
  {
    return std::log(x/(1-x));
  }

  template <typename T>
  T ilogit(const T& x)
  {
    return std::exp(x)/(1+std::exp(x));
  }
}


#endif
