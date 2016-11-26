#ifndef PTSAMPLE_HPP
#define PTSAMPLE_HPP

#include "ensemble_sample.hpp"

namespace mcmc_utilities
{
  template <typename T_logprob,typename T_var>
  typename std::result_of<T_logprob(T_var)>::type exchange_prob(T_logprob&& logprob,
								const T_var& var1,
								const T_var& var2,
								typename std::result_of<T_logprob(T_var)>::type beta1,
								typename std::result_of<T_logprob(T_var)>::type beta2)
  {
    using T=typename std::result_of<T_logprob(T_var)>::type;
    //beta2<beta1 2 is a hi-T system
    //if logprob(var2)<logprob(var1) => logprob(var2)-logprob(var1)<0
    //exp(...)<1
    T p1=logprob(var1);
    T p2=logprob(var2);
    if(std::isinf(p1)||std::isinf(p2))
      {
	//return 0;
      }
    return std::min(as<T>(1),std::exp((beta2-beta1)*(-p2+p1)));
  }

  template <typename T_container>
  void swap(T_container& arr,size_t i1,size_t i2)
  {
    auto temp=clone<typename element_type_trait<T_container>::element_type>(get_element(arr,i1));
    set_element(arr,i1,get_element(arr,i2));
    set_element(arr,i2,temp);
  }

  template <typename T_container,typename T_rng>
  void shuffle(T_container& arr,
	       T_rng& rng)
  {
    size_t n=get_size(arr);
    //for(int i=n-1;i>=0;--i)
    for(size_t j=0;j<n;++j)
      {
	int i=n-1-j;
	int i2=0;
	do
	  {
	    i2=rng()*n;
	  }
	while(i2==n);
	//std::swap(get_element(arr,i),get_element(arr,i2));
	swap(arr,i,i2);
      }
  }

  template <typename T_logprob,typename T_ensemble_list,typename T_beta_list,typename T_rng>
  T_ensemble_list ptsample(T_logprob&& logprob,
			   const T_ensemble_list& ensemble_list,
			   T_rng&& rng,
			   const T_beta_list& beta_list,
			   bool perform_swap,
			   size_t nthread_allowed=1,
			   typename std::result_of<T_logprob(typename element_type_trait<typename element_type_trait<T_ensemble_list>::element_type>::element_type)>::type a=2)
  {
    using T=typename std::result_of<T_logprob(typename element_type_trait<typename element_type_trait<T_ensemble_list>::element_type>::element_type)>::type;
    using T_var=typename element_type_trait<typename element_type_trait<T_ensemble_list>::element_type>::element_type;
    auto new_ensemble_list=clone(ensemble_list);
    size_t ntemp=get_size(ensemble_list);
    size_t nwalker=get_size(get_element(ensemble_list,0));
    
    if(perform_swap)
      {
	/*
	for(size_t i=0;i<ntemp;++i)
	  {
	    shuffle(get_element(new_ensemble_list,i),rng);
	  }
	*/
	for(size_t i=0;i<ntemp-1;++i)
	  {
	    T beta1=beta_list[i];
	    T beta2=beta_list[i+1];

	    if(beta1==beta2)
	      {
		mcmc_exception e("beta list should not contain duplicated elements");
		throw e;
	      }
	    
	    for(size_t j=0;j<nwalker;++j)
	      {
		auto var1=as<T_var>(get_element(get_element(new_ensemble_list,i),j));
		auto var2=as<T_var>(get_element(get_element(new_ensemble_list,i+1),j));
		T ep=exchange_prob(logprob,var1,var2,beta1,beta2);
		if(urng<T>(rng)<ep)
		  {
		    auto temp=clone<T_var>(get_element(get_element(new_ensemble_list,i),j));
		    set_element(get_element(new_ensemble_list,i),j,
				get_element(get_element(new_ensemble_list,i+1),j));
		    set_element(get_element(new_ensemble_list,i+1),j,temp);
		  }
	      }
	  }
      }
    
    
    for(size_t i=0;i<ntemp;++i)
      {
	T beta=get_element(beta_list,i);
	set_element(new_ensemble_list,i,ensemble_sample([&logprob,beta](const T_var& x){
	      T lp=logprob(x);
	      if(std::isinf(lp))
		{
		  return lp;
		}
	      return lp*beta;
	    },get_element(new_ensemble_list,i),rng,nthread_allowed,a));
      }
    return new_ensemble_list;
  }


  template <typename T_prior,typename T_likelihood,typename T_ensemble_list,typename T_beta_list,typename T_rng>
  T_ensemble_list ptsample_pl(T_prior&& logprior,
			   T_likelihood&& loglikelihood,
			   const T_ensemble_list& ensemble_list,
			   T_rng&& rng,
			   const T_beta_list& beta_list,
			   bool perform_swap,
			   size_t nthread_allowed=1,
			   typename std::result_of<T_prior(typename element_type_trait<typename element_type_trait<T_ensemble_list>::element_type>::element_type)>::type a=2)
  {
    using T=typename std::result_of<T_prior(typename element_type_trait<typename element_type_trait<T_ensemble_list>::element_type>::element_type)>::type;
    using T_var=typename element_type_trait<typename element_type_trait<T_ensemble_list>::element_type>::element_type;
    auto new_ensemble_list=clone(ensemble_list);
    size_t ntemp=get_size(ensemble_list);
    size_t nwalker=get_size(get_element(ensemble_list,0));
    
    if(perform_swap)
      {
	/*
	for(size_t i=0;i<ntemp;++i)
	  {
	    shuffle(get_element(new_ensemble_list,i),rng);
	  }
	*/
	for(size_t i=0;i<ntemp-1;++i)
	  {
	    T beta1=beta_list[i];
	    T beta2=beta_list[i+1];

	    if(beta1==beta2)
	      {
		mcmc_exception e("beta list should not contain duplicated elements");
		throw e;
	      }
	    
	    for(size_t j=0;j<nwalker;++j)
	      {
		auto var1=as<T_var>(get_element(get_element(new_ensemble_list,i),j));
		auto var2=as<T_var>(get_element(get_element(new_ensemble_list,i+1),j));
		T ep=exchange_prob([&logprior,&loglikelihood](const T_var& x){return logprior(x)+loglikelihood(x);},var1,var2,beta1,beta2);
		if(urng<T>(rng)<ep)
		  {
		    auto temp=clone<T_var>(get_element(get_element(new_ensemble_list,i),j));
		    set_element(get_element(new_ensemble_list,i),j,
				get_element(get_element(new_ensemble_list,i+1),j));
		    set_element(get_element(new_ensemble_list,i+1),j,temp);
		  }
	      }
	  }
      }
    
    
    for(size_t i=0;i<ntemp;++i)
      {
	T beta=get_element(beta_list,i);
	set_element(new_ensemble_list,i,ensemble_sample([&logprior,&loglikelihood,beta](const T_var& x){
	      T lp=logprior(x);
	      T ll=loglikelihood(x);
	      if(std::isinf(lp)||std::isinf(ll))
		{
		  return lp;
		}
	      return lp+ll*beta;
	    },get_element(new_ensemble_list,i),rng,nthread_allowed,a));
      }
    return new_ensemble_list;
  }
}

#endif
