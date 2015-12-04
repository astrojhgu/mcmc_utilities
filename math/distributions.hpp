#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP
#include <cmath>
#include <cassert>
#include "../core/error_handler.hpp"
#include "functions.hpp"
namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  T_p logdgamma(T_var x,T_p r,T_p lambda)
  {
    if(x<0||lambda<=0||r<=0)
      {
	throw var_out_of_range();
      }
    return r*std::log(lambda)+(r-1)*std::log(x)-lambda*x-lgamma(r);
  }

  template <typename T_p,typename T_var>
  T_p logdt(T_var x,T_p mu,T_p tau,T_p k)
  {
    if(tau<=0||k<=0)
      {
	throw var_out_of_range();
      }
    static const T_p pi=std::atan(1)*4;
    return lgamma((k+1)/2)-lgamma(k/2)+std::log(tau/k/pi)/2-(k+1)/2*std::log(1+tau*(x-mu)*(x-mu)/k);
  }

  template <typename T_p,typename T_var>
  T_p logdnorm(T_var x,T_p mu,T_p sigma)
  {
    if(sigma<=0)
      {
	std::cerr<<"sigma="<<sigma<<std::endl;
	std::cerr<<"x="<<x<<std::endl;
	//assert(0);
	throw var_out_of_range();
      }
    static const T_p pi=std::atan(1)*4;
    //return std::log(tau/(2*pi))/2-tau*(x-mu)*(x-mu)/2;
    return -std::log(2*pi*sigma*sigma)/2-(x-mu)*(x-mu)/(2*sigma*sigma);
  }

  template <typename T_p,typename T_var>
  T_p logdlnorm(T_var x,T_p mu,T_p sigma)
  {
    if(sigma<=0||x<=0)
      {
	throw var_out_of_range();
      }
    static const T_p pi=std::atan(1)*4;
    T_p result=-std::log(2*pi*sigma*sigma)/2-std::log(x)-(std::log(x)-mu)*(std::log(x)-mu)/(2.*sigma*sigma);
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
    T_p result=std::log(alpha)+alpha*std::log(c)-(alpha+1)*std::log(x);
    return result;
  }


  template <typename T_p,typename T_var>
  T_p logdbivnorm(const T_var& x,
		  const T_p& mu1,const T_p& mu2,
		  const T_p& sigma11,const T_p& sigma22,T_p& sigma12)
		  
  {
    static const T_p pi=std::atan(1)*4;
    
    T_p sigma21=sigma12;
    
    T_p sigma1_sq=sigma11*sigma11+sigma12*sigma12;
    T_p sigma2_sq=sigma21*sigma21+sigma22*sigma22;

    T_p sigma1=std::sqrt(sigma1_sq);
    T_p sigma2=std::sqrt(sigma2_sq);

    T_p rho=(sigma11*sigma21+sigma12*sigma22)/(sigma1*sigma2);
      
    T_p x1=x[0];
    T_p x2=x[1];
    T_p z=(x1-mu1)*(x1-mu1)/sigma1_sq+
      (x2-mu2)*(x2-mu2)/sigma2_sq-
      2*rho*(x1-mu1)*(x2-mu2)/(sigma1*sigma2);
    return (-z/(2*(1-rho*rho)))-std::log(2*pi*sigma1*sigma2*std::sqrt(1-rho*rho));

  }

}


#endif
