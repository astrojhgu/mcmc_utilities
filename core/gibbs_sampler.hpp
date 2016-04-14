#ifndef GIBBS_SAMPLER_HPP
#define GIBBS_SAMPLER_HPP
//#include "rejection_sampler_1d.hpp"
#include "base_urand.hpp"
//#include "arms.hpp"
#include "discrete_sample.hpp"
#include "continuous_sample.hpp"
#include <vector>
#include <memory>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  void gibbs_sample(const probability_density_md<T_p,T_var>& pd,T_var& init_var,base_urand<T_p>& rnd)
  {
    size_t idx=0;
    typedef typename element_type_trait<T_var>::element_type T_var1;
    
    T_var1 xprev=0.;
    for(idx=0;idx<get_size(init_var);++idx)
      {
	xprev=get_element(init_var,idx);
		
	//T_var1 x=sampler.sample(cpd,xprev);
	//T_var1 x=arms(cpd,xprev,10,rnd);
	auto var_range=pd.var_range(init_var,idx);

	std::vector<T_var1> xinit(pd.init_points(init_var,idx));
	if(xinit.size()==0)
	  {
	    xinit.resize(5);
	    T_var1 xl=var_range.first;
	    T_var1 xr=var_range.second;

	    for(size_t n1=0;n1<xinit.size();++n1)
	      {
		xinit[n1]= xl+(xr-xl)/(xinit.size()+1)*(n1+1);
	      }
	  }
	
	T_var1 x=continuous_sample([&](const T_var1& x){
	    set_element(init_var,idx,x);
	    T_p result= pd.eval_log(init_var,idx);
	    if(!std::isfinite(result))
	      {
		std::cerr<<"inside gibbs sampler:"<<std::endl;
		std::cerr<<"x="<<x<<std::endl;
		std::cerr<<"init_var["<<idx<<"]="<<init_var[idx]<<std::endl;
	      }
	    return result;
	  },var_range,xinit,xprev,10,rnd);
	set_element(init_var,idx,x);
      }
  }

  template <typename T_p,typename T_var,typename T_urand>
  void gibbs_sample1(const probability_density_md<T_p,T_var>& pd,
		    T_var& init_var,size_t idx,base_urand<T_p>& rnd)
  {
    for(int i=0;i<get_size(init_var);++i)
      {
	auto xrange=pd.var_range(init_var,i);
	if(xrange.first>get_element(init_var,i)||xrange.second<get_element(init_var,i))
	  {
	    var_out_of_range e;
	    e.attach_message("gibbs sampler #1");
	    throw e;
	  }
      }
    
    if(idx>=get_size(init_var))
      {
	throw index_out_of_range();
      }
    typedef typename element_type_trait<T_var>::element_type T_var1;
    
    T_var1 xprev=0.;


    xprev=get_element(init_var,idx);
    auto var_range=pd.var_range(init_var,idx);    
    std::vector<T_var1> xinit(pd.init_points(init_var,idx));
    if(xinit.size()==0)
      {
	xinit.resize(5);
	
	
	T_var1 xl=var_range.first;
	T_var1 xr=var_range.second;
	
	for(size_t n1=0;n1<xinit.size();++n1)
	  {
	    xinit[n1]= xl+(xr-xl)/(xinit.size()+1)*(n1+1);
	  }
      }

    T_var1 x=continuous_sample([&](const T_var1& x){
	set_element(init_var,idx,x);
	T_p result= pd.eval_log(init_var,idx);
	return result;
      }, var_range,xinit,xprev,10,rnd);
    set_element(init_var,idx,x);
  }
  
};

#endif
