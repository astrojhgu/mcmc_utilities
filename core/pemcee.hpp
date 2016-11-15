#ifndef PEMCEE_HPP
#define PEMCEE_HPP
#include <cassert>
#include <cmath>
#include <functional>
#include <type_traits>
#include <thread>
#include <mutex>
#include <atomic>
#include <sstream>
#include "mcmc_traits.hpp"
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
  
  template <typename T_logprob,typename T_ensemble>
  const T_ensemble pemcee(const T_logprob& logprob,
	  const T_ensemble& ensemble,
	  base_urand<typename std::result_of<T_logprob(typename T_ensemble::value_type)>::type>& rnd,
	  typename std::result_of<T_logprob(typename T_ensemble::value_type)>::type a=2)
  {
    const size_t K=ensemble.size();
    const size_t n=get_element(ensemble,0).size();
    const size_t half_K=K/2;
    using T=typename std::result_of<T_logprob(typename T_ensemble::value_type)>::type;
    using T_var=typename T_ensemble::value_type;
    if(K%2!=0)
      {
	throw mcmc_exception("number of walkers must be even");
      }
    T_ensemble ensemble_half(ensemble);

    auto task=[&](size_t k)
      {
	const size_t i=k<half_K?0:1;
	const size_t ni=1-i;
	size_t j=0;
	do
	  {
	    j=rnd()*(half_K);
	  }
	while(j==half_K);
	T z=draw_z(rnd,a);
	T_var Y(get_element(ensemble,k));
	for(int l=0;l<Y.size();++l)
	  {
	    T y=get_element(get_element(ensemble,j+half_K*ni),l)+z*(get_element(get_element(ensemble,k),l)-get_element(get_element(ensemble,j+half_K*ni),l));
	    set_element(Y,l,y);
	    if(std::isnan(y)||std::isinf(y))
	      {
		nan_or_inf e;
		e.attach_message("inf or nan for y\n");
		std::ostringstream oss;
		oss<<"y="<<y<<std::endl;
		e.attach_message(oss.str());
	      }
	  }
	T q=std::exp((n-1)*std::log(z)+(logprob(Y)-logprob(get_element(ensemble,k))));
	if(std::isnan(q)||std::isinf(q))
	  {
	    nan_or_inf e;
	    e.attach_message("inf or nan\n");
	    
	    std::ostringstream oss;
	    oss<<"q="<<q<<std::endl;
	    oss<<"Y=";
	    for(int l=0;l<n;++l)
	      {
		oss<<get_element(Y,l)<<" ";
	      }
	    oss<<"\nlogprob(Y)="<<logprob(Y)<<"\n";

	    oss<<"X=";
	    for(int l=0;l<n;++l)
	      {
		oss<<get_element(get_element(ensemble,k),l)<<" ";
	      }
	    oss<<"\nlogprob(X)="<<logprob(get_element(ensemble,k))<<"\n";
	    e.attach_message(oss.str());
	    throw e;
	  }
	
	T r=rnd();
	if(r<=q)
	  {
	    set_element(ensemble_half,k,Y);
	  }
      };

    if(!rnd.is_parallel())
      {
	for(size_t k=0;k<K;++k)
	  {
	    task(k);
	  }
      }
    else
      {
	std::vector<std::thread> pool;
	pool.reserve(K);
	size_t nthread_allowed=std::thread::hardware_concurrency();
	size_t nrunning_thread(0);
	size_t k(0);
	std::mutex mx;
	std::function<void()> chain_task([&nrunning_thread,&k,&task,nthread_allowed,&pool,&chain_task,K,&mx](){
	    mx.lock();
	    int k1=k;
	    //std::cerr<<"k="<<k<<std::endl;
	    nrunning_thread++;
	    k++;
	    mx.unlock();
	    task(k1);
	    mx.lock();
	    if(nrunning_thread<nthread_allowed&&pool.size()<K)
	      {
		pool.push_back(std::thread(chain_task));
		//assert(pool.size()<=K);
	      }
	    nrunning_thread--;
	    mx.unlock();
	  });
	mx.lock();
	for(int i=0;i<nthread_allowed;++i)
	  {
	    pool.push_back(std::thread(chain_task));
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
