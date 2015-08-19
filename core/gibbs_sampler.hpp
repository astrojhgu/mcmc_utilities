#ifndef GIBBS_SAMPLER_HPP
#define GIBBS_SAMPLER_HPP
//#include "rejection_sampler_1d.hpp"
#include "uvsampler.hpp"
#include "discrete_sample.hpp"
#include <vector>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  void gibbs_sample(const probability_density_md<T_p,T_var>& pd,T_var& init_var,const uvsampler<T_p,typename element_type_trait<T_var>::element_type>& sampler)
  {
    size_t idx=0;
    typedef typename element_type_trait<T_var>::element_type T_var1;
    class conditional_probability_density
      :public probability_density_1d<T_p,typename element_type_trait<T_var>::element_type>
    {
    public:
      size_t* p_idx;
      T_var* p_init_var;
      const probability_density_md<T_p,T_var>* ppd;
    private:
      T_p do_eval_log(const T_var1& x1)const
      {
	set_element(*p_init_var,*p_idx,x1);
	T_p result= ppd->eval_log(*p_init_var);
	return result;
      }
      
      std::pair<T_var1,T_var1> do_var_range()const
      {
	//ppd->var_range(xmin,xmax,*p_init_var,(*p_idx));
	return ppd->var_range(*p_init_var,*p_idx);
      }

      std::vector<T_var1> do_init_points()const
      {
	//int n=ppd->num_init_points(*p_init_var,*p_idx);
	std::vector<T_var1> xinit(ppd->init_points(*p_init_var,*p_idx));
	if(xinit.size()==0)
	  {
	    xinit.resize(3);
	    
	    std::pair<T_var1,T_var1> xrange(this->var_range());
	    T_var1 xl=xrange.first;
	    T_var1 xr=xrange.second;

	    for(int n1=0;n1<3;++n1)
	      {
		if(n1!=1)
		  {
		    xinit[n1]= xl+(xr-xl)/(3+1)*(n1+1);
		  }
		else
		  {
		    xinit[n1]= find_peak(*this);
		  }
	      }
	  }
	return xinit;
      }

      std::vector<T_var1> do_candidate_points()const
      {
	return ppd->candidate_points(*p_init_var,*p_idx);
      }
    }cpd;
    cpd.p_idx=&idx;
    cpd.p_init_var=&init_var;
    cpd.ppd=&pd;
    T_var1 xprev=0.;
    for(idx=0;idx<get_size(init_var);++idx)
      {
	xprev=get_element(init_var,idx);
	T_var1 x=sampler.sample(cpd,xprev);
	set_element(init_var,idx,x);
      }
  }

  template <typename T_p,typename T_var,typename T_urand>
  void gibbs_sample1(const probability_density_md<T_p,T_var>& pd,
		    T_var& init_var,size_t idx,const uvsampler<T_p,T_var>& sampler)
  {
    if(idx>=get_size(init_var))
      {
	throw index_out_of_range();
      }
    typedef typename element_type_trait<T_var>::element_type T_var1;
    class conditional_probability_density
      :public probability_density_1d<T_p,typename element_type_trait<T_var>::element_type>
    {
    public:
      size_t* p_idx;
      T_var* p_init_var;
      const probability_density_md<T_p,T_var>* ppd;
    private:
      T_p do_eval_log(const T_var1& x1)const
      {
	set_element(*p_init_var,*p_idx,x1);
	T_p result= ppd->eval_log(*p_init_var);
	return result;
      }
 
      std::pair<T_var1,T_var1> do_var_range()const
      {
	//ppd->var_range(xmin,xmax,*p_init_var,(*p_idx));
	return ppd->var_range(*p_init_var,*p_idx);
      }

      std::vector<T_var1> do_init_points()const
      {
	//int n=ppd->num_init_points(*p_init_var,*p_idx);
	std::vector<T_var1> xinit(ppd->init_points(*p_init_var,*p_idx));
	if(xinit.size()==0)
	  {
	    xinit.resize(3);
	    
	    std::pair<T_var1,T_var1> xrange(this->var_range());
	    T_var1 xl=xrange.first;
	    T_var1 xr=xrange.second;

	    for(int n1=0;n1<3;++n1)
	      {
		if(n1!=1)
		  {
		    xinit[n1]= xl+(xr-xl)/(3+1)*(n1+1);
		  }
		else
		  {
		    xinit[n1]= find_peak(*this);
		  }
	      }
	  }
	return xinit;
      }

      std::vector<T_var1> do_candidate_points()const
      {
	return ppd->candidate_points(*p_init_var,*p_idx);
      }
      
    }cpd;
    cpd.p_idx=&idx;
    cpd.p_init_var=&init_var;
    cpd.ppd=&pd;
    T_var1 xprev=0.;


    xprev=get_element(init_var,idx);
    T_var1 x=sampler.sample(cpd,xprev);
    set_element(init_var,idx,x);
  }
  
};

#endif
