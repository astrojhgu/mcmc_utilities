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

  template <typename T>
  T logdgamma(T x,T r,T lambda)
  {
    if(x<0||lambda<=0||r<=0)
      {
	throw var_out_of_range();
      }
    return r*std::log(lambda)+(r-1)*std::log(x)-lambda*x-lgamma(r);
  }

  template <typename T>
  T logdt(T x,T mu,T tau,T k)
  {
    if(tau<=0||k<=0)
      {
	throw var_out_of_range();
      }
    static const T pi=std::atan(1)*4;
    return lgamma((k+1)/2)-lgamma(k/2)+std::log(tau/k/pi)/2-(k+1)/2*std::log(1+tau*(x-mu)*(x-mu)/k);
  }

  template <typename T>
  T logdnorm(T x,T mu,T tau)
  {
    if(tau<=0)
      {
	throw var_out_of_range();
      }
    static const T pi=std::atan(1)*4;
    return std::log(tau/(2*pi))/2-tau*(x-mu)*(x-mu)/2;
  }

  template <typename T_p,typename T_var>
  T_p log_factorial(const T_var& n)
  {
    if(n==0)
      {
	return 0;
      }
    T_p result=0;
    
    for(T_var x=1;x<=n;++x)
      {
	result+=std::log(x);
      }
    return result;
  }

  template <typename T_p,typename T_var>
  T_p log_CN(const T_var& m,const T_var& n)
  {
    T_p result= log_factorial<T_p,T_var>(m)-log_factorial<T_p,T_var>(n)-log_factorial<T_p,T_var>(m-n);
    return result;
  }

  template <typename T_p,typename T_var>
  T_p logdbin(const T_var& x,const T_p& p,const T_var& n)
  {
    T_p result= log_CN<T_p,T_var>(n,x)+x*std::log(p)+(n-x)*std::log(1-p);
    return result;
  }
}


#endif
