#ifndef PEMCEE_HPP
#define PEMCEE_HPP
#include <cassert>
#include <cmath>
#include "base_urand.hpp"
#include "distribution.hpp"
#include "error_handler.hpp"

/*
  implemeted the emcee algorithm according to
  http://adsabs.harvard.edu/abs/2013PASP..125..306F
  https://arxiv.org/abs/1202.3665
*/

namespace mcmc_utilities
{
  template <typename T>
  T draw_z(base_urand<T>& rnd,T a)
  {
    const T sqrt_a=std::sqrt(a);
    T p=rnd()*2*(sqrt_a-1/sqrt_a);
    return std::pow(static_cast<T>(1/sqrt_a+p/2),static_cast<T>(2));
  }
  
  template <typename T,template <typename TE> class T_vector,typename T_var>
  void pemcee(const probability_density_md<T,T_var,T_vector>& prob,
	      T_vector<T_var>& ensemble,
	      base_urand<T>& rnd,
	      T a=2)
  {
    size_t K=ensemble.size();
    size_t n=ensemble[0].size();
    if(K%2!=0)
      {
	throw mcmc_exception("number of positions must be even");
      }
    T_vector<T_var> ensemble_half(ensemble);
    for(size_t i:{0,1})
      {
	size_t ni=1-i;
	for(size_t k=0;k<K/2;++k)
	  {
	    size_t j=0;
	    do
	      {
		j=rnd()*(K/2);
	      }
	    while(j==K/2);
	    T z=draw_z(rnd,a);
	    T_var Y(ensemble[k+K/2*i]);
	    for(int l=0;l<Y.size();++l)
	      {
		Y[l]=ensemble[j+K/2*ni][l]+z*(ensemble[k+K/2*i][l]-ensemble[j+K/2*ni][l]);
	      }
	    T q=std::exp((n-1)*std::log(z)+(prob.eval_log(Y)-prob.eval_log(ensemble[k+K/2*i])));
	    T r=rnd();
	    if(r<=q)
	      {
		ensemble_half[k+K/2*i]=Y;
	      }
	  }
      }
    ensemble.swap(ensemble_half);
  }
}

#endif
