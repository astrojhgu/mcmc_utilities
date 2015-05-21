#ifndef MCMC_BAS_UTIL
#define MCMC_BAS_UTIL
#define OPT_HEADER
#include "mcmc_traits.hpp"
#include <algorithm>
#include <limits>
#include "distribution.hpp"
namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  struct dist_adapter
  {
      T_var xl,xr;
      const probability_density_1d<T_p,T_var>* ppd;

      T_p operator()(const T_var& x)const
      {
	if(x<=xl||x>=xr)
	  {
	    return std::numeric_limits<T_p>::max();
	  }
	return -(ppd->eval_log(x));
      }
  };
  
  template <typename T>
  T sqr(T x)
  {
    return x*x;
  }
  
  
  template <typename T>
  void shft3(T&a,T& b,T& c,T d)
  {
    mcmc_assign(a,b);
    mcmc_assign(b,c);
    mcmc_assign(c,d);
  }

  template <typename T>
  void shft(T& a,T& b,T& c,T d)
  {
    mcmc_assign(a,b);
    mcmc_assign(b,c);
    mcmc_assign(c,d);
  }
  //  template <typename T>  
  //  void swap(T& ax,T& bx)
  //{
    //  swap(ax,bx);
    //    T temp;
    //mcmc_assign(temp,ax);
    //mcmc_assign(ax,bx);
    //mcmc_assign(bx,temp);
    //}
  
  template <typename T>
  T sign(const T& a,const T& b)
  {
    return b>=0?T(a>=0?T(a):T(-a)):T(a>=0?T(-a):T(a));
  }


  template <typename T>
  void mov3(T& a,T& b,T& c, T& d,T& e,T& f)
  {
    mcmc_assign(a,d);mcmc_assign(b,e);mcmc_assign(c,f);
  }
}

#endif