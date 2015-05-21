#ifndef SPECIAL_FUNCTION_HPP
#define SPECIAL_FUNCTION_HPP
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
  T logdnorm(T x,T mu,T sigma)
  {
    if(sigma<=0)
      {
	assert(0);
	throw var_out_of_range();
      }
    static const T pi=std::atan(1)*4;
    //return std::log(tau/(2*pi))/2-tau*(x-mu)*(x-mu)/2;
    return -std::log(2*pi*sigma*sigma)/2-(x-mu)*(x-mu)/(2*sigma*sigma);
  }

  template <typename T>
  T logdlnorm(T x,T mu,T sigma)
  {
    if(sigma<=0||x<=0)
      {
	assert(0);
	throw var_out_of_range();
      }
    static const T pi=std::atan(1)*4;
    T result=-std::log(2*pi*sigma*sigma)/2-std::log(x)-(std::log(x)-mu)*(std::log(x)-mu)/(2.*sigma*sigma);
    assert(!std::isnan(result)&&!std::isinf(result));
    return result;
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

  template <typename T_p,typename T_var>
  T_p logdpoisson(const T_var& x,const T_p& lambda)
  {
    T_p result=x*std::log(lambda)-lambda-log_factorial<T_p,T_var>(x);
    return result;
  }

  template <typename T_p,typename T_var>
  T_p logdpar(const T_var& x,const T_var& c,const T_p& alpha)
  {
    if(x<c||c<=0)
      {
	throw var_out_of_range();
      }
    assert(x>0);
    T_p result=std::log(alpha)+alpha*std::log(c)-(alpha+1)*std::log(x);
    assert(!std::isnan(result)&&!std::isinf(result));
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
