#ifndef ENSEMBLE_SAMPLER_HPP
#define ENSEMBLE_SAMPLER_HPP
#include <cassert>
#include <cmath>
#include <functional>
#include <type_traits>
#include <thread>
#include <mutex>
#include <atomic>
#include <sstream>
#include "mcmc_traits.hpp"
#include "error_handler.hpp"

/*
  implemeted the emcee algorithm according to
  http://adsabs.harvard.edu/abs/2013PASP..125..306F
  https://arxiv.org/abs/1202.3665
*/

namespace mcmc_utilities
{
  template <typename T,typename T_rng>
  T draw_z(T_rng&& rng,T a)
  {
    const T sqrt_a=std::sqrt(a);
    T p=urng<T>(rng)*2*(sqrt_a-1/sqrt_a);
    return std::pow(static_cast<T>(1/sqrt_a+p/2),static_cast<T>(2));
  }
  
  template <typename T_logprob,typename T_ensemble,typename T_rng>
  T_ensemble ensemble_sample(T_logprob&& logprob,
			     const T_ensemble& ensemble,
			     T_rng&& rng,
			     size_t nthread_allowed=1,
			     typename std::result_of<T_logprob(typename element_type_trait<T_ensemble>::element_type)>::type a=2)
  {
    using T=typename std::result_of<T_logprob(typename element_type_trait<T_ensemble>::element_type)>::type;
    using T_var=typename element_type_trait<T_ensemble>::element_type;
    
    const size_t K=get_size(ensemble);
    
    if(K==0)
      {
	throw mcmc_exception("number of walkers must not be zero");
      }
    if(K%2!=0)
      {
	throw mcmc_exception("number of walkers must be even");
      }

    std::vector<std::vector<size_t> > walker_group{std::vector<size_t>(),std::vector<size_t>()};
    std::vector<size_t> walker_group_id(K);
    for(auto& i:walker_group)
      {
	i.reserve(K/2);
      }
    for(size_t i=0;i<K;++i)
      {
	size_t gid=urng<T>(rng)<static_cast<T>(0.5)?0:1;
	if(walker_group[gid].size()==K/2)
	  {
	    gid=1-gid;
	  }
	walker_group[gid].push_back(i);
	walker_group_id[i]=gid;
      }

    
    nthread_allowed=nthread_allowed<1?1:nthread_allowed;
    const size_t n=get_size(get_element(ensemble,0));
    const size_t half_K=K/2;
    //T_ensemble ensemble_half(ensemble);
    auto ensemble_half(clone(ensemble));

    auto task=[&](size_t k,typename std::remove_reference<T_logprob>::type& logprob)
      {
	const size_t i=walker_group_id[k];
	const size_t ni=1-i;
	size_t j=0;
	do
	  {
	    j=urng<T>(rng)*(half_K);
	  }
	while(j>=half_K);
	T z=draw_z(rng,a);
	//T_var Y(get_element(ensemble,k));
	T_var Y(clone(get_element(ensemble,k)));
	for(size_t l=0;l<get_size(Y);++l)
	  {
	    T y=as<T>(get_element(get_element(ensemble,walker_group[ni][j]),l))+z*(as<T>(get_element(get_element(ensemble,k),l))-as<T>(get_element(get_element(ensemble,walker_group[ni][j]),l)));
	    if(std::isnan(y)||std::isinf(y))
	      {
		nan_or_inf e;
		e.attach_message("inf or nan for y\n");
		std::ostringstream oss;
		oss<<"y="<<y<<std::endl;
		e.attach_message(oss.str());
		
	      }
	    set_element(Y,l,as<typename element_type_trait<T_var>::element_type>(y));
	  }
	T lpY=logprob(Y);
	T lpLastY=logprob(get_element(ensemble,k));
	if(std::isinf(lpLastY))
	  {
	    nan_or_inf e;
	    e.attach_message("inf or nan\n");
	    e.attach_message("the logprob of the members in last ensemble should not be inf");
	    std::ostringstream oss;
	    oss<<"last Y=";
	    for(size_t l=0;l<n;++l)
	      {
		oss<<as<T>(get_element(get_element(ensemble,k),l))<<" ";
	      }
	    oss<<"\nlogprob(Y)="<<lpLastY<<"\n";
	    e.attach_message(oss.str());
	    throw e;
	  }

	if(!std::isinf(lpY)&&!std::isnan(lpY))
	  {
	    T q=std::exp((n-1)*std::log(z)+lpY-lpLastY);
	    if(std::isnan(q))
	      {
		nan_or_inf e;
		e.attach_message("inf or nan\n");
		
		std::ostringstream oss;
		oss<<"q="<<q<<std::endl;
		oss<<"Y=";
		for(size_t l=0;l<n;++l)
		  {
		    oss<<as<T>(get_element(Y,l))<<" ";
		  }
		oss<<"\nlogprob(Y)="<<logprob(Y)<<"\n";
		
		oss<<"X=";
		for(size_t l=0;l<n;++l)
		  {
		    oss<<as<T>(get_element(get_element(ensemble,k),l))<<" ";
		  }
		oss<<"\nlogprob(X)="<<logprob(get_element(ensemble,k))<<"\n";
		oss<<"z="<<z<<std::endl;
		e.attach_message(oss.str());
		
		throw e;
	      }
	    
	    T r=urng<T>(rng);
	    if(r<=q)
	      {
		set_element(ensemble_half,k,as<typename element_type_trait<T_ensemble>::element_type>(Y));
	      }
	  }
      };

    if(nthread_allowed==1)
      {
	for(size_t k=0;k<K;++k)
	  {
	    task(k,logprob);
	  }
      }
    else
      {
	std::vector<typename std::remove_reference<T_logprob>::type> logprob_vec;
	logprob_vec.reserve(nthread_allowed);
	for(size_t i=0;i<nthread_allowed;++i)
	  {
	    logprob_vec.push_back(clone(logprob));
	  }
	
	std::vector<std::thread> pool;
	pool.reserve(K);
	size_t nrunning_thread(0);
	size_t k(0);
	std::mutex mx;
	std::function<void(typename std::remove_reference<T_logprob>::type&)> chain_task([&nrunning_thread,&k,&task,nthread_allowed,&pool,&chain_task,K,&mx](typename std::remove_reference<T_logprob>::type& logprob){
	    mx.lock();
	    size_t k1=k;
	    //std::cerr<<"k="<<k<<std::endl;
	    nrunning_thread++;
	    k++;
	    mx.unlock();
	    task(k1,logprob);
	    mx.lock();
	    if(nrunning_thread<nthread_allowed&&pool.size()<K)
	      {
		pool.push_back(std::thread(chain_task,std::ref(logprob)));
		//assert(pool.size()<=K);
	      }
	    nrunning_thread--;
	    mx.unlock();
	  });
	mx.lock();
	for(size_t i=0;i<nthread_allowed;++i)
	  {
	    pool.push_back(std::thread(chain_task,std::ref(logprob_vec[i])));
	  }
	mx.unlock();
	//std::cerr<<pool.size()<<std::endl;
	bool any_running=false;
	do
	  {
	    any_running=false;
	    for(auto& t:pool)
	      {
		if(t.joinable())
		  {
		    any_running=true;
		    t.join();
		  }
	      }
	  }
	while(any_running||pool.size()<K);
      }
    
    return (ensemble_half);
  }
}

#endif
