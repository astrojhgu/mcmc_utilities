#ifndef GIBBS_SAMPLER_HPP
#define GIBBS_SAMPLER_HPP
//#include "rejection_sampler_1d.hpp"
#include "arms.hpp"
#include <vector>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var,typename T_urand>
  void gibbs_sample(const probability_density_md<T_p,T_var>& pd,T_var& init_var,bool dometrop,const T_urand& urand, int sample_cnt=10)
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
      
      void do_var_range(T_var1& xmin,T_var1& xmax)const
      {
	ppd->var_range(xmin,xmax,*p_init_var,(*p_idx));
	T_var1 x=get_element(*p_init_var,*p_idx);
      }

      size_t do_num_init_points()const
      {
	int n=ppd->num_init_points(*p_init_var,*p_idx);
	if(n==0)
	  {
	    return 3;
	  }
	else
	  {
	    return n;
	  }
      }

      T_var1 do_init_point(size_t n1)const
      {
	int n=ppd->num_init_points(*p_init_var,*p_idx);
	if(n==0)
	  {
	    T_var1 xl,xr;
	    this->var_range(xl,xr);
	    if(n1!=1)
	      {
		return xl+(xr-xl)/(this->num_init_points()+1)*(n1+1);
	      }
	    return find_peak(*this);
	  }
	else
	  {
	    return ppd->init_point(n1,*p_init_var,*p_idx);
	  }
      }
    }cpd;
    cpd.p_idx=&idx;
    cpd.p_init_var=&init_var;
    cpd.ppd=&pd;
    T_var1 xprev=0.;
    std::vector<T_var1> xsamp(sample_cnt);
    for(idx=0;idx<get_size(init_var);++idx)
      {
	xprev=get_element(init_var,idx);
	arms_simple(cpd,xprev,xsamp,dometrop,urand);
	set_element(init_var,idx,xsamp.back());
      }
  }

  template <typename T_p,typename T_var,typename T_urand>
  void gibbs_sample(const probability_density_md<T_p,T_var>& pd,
		    T_var& init_var,size_t idx,bool dometrop,const T_urand& urand, int sample_cnt=10)
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
      
      void do_var_range(T_var1& xmin,T_var1& xmax)const
      {
	ppd->var_range(xmin,xmax,*p_init_var,(*p_idx));
	T_var1 x=get_element(*p_init_var,*p_idx);
      }


      size_t do_num_init_points()const
      {
	int n=ppd->num_init_points(*p_init_var,*p_idx);
	if(n==0)
	  {
	    return 3;
	  }
	else
	  {
	    return n;
	  }
      }

      T_var1 do_init_point(size_t n1)const
      {
	int n=ppd->num_init_points(*p_init_var,*p_idx);
	if(n==0)
	  {
	    T_var1 xl,xr;
	    this->var_range(xl,xr);
	    if(n1!=1)
	      {
		return xl+(xr-xl)/(this->num_init_points()+1)*(n1+1);
	      }
	    return find_peak(*this);
	  }
	else
	  {
	    return ppd->init_point(n1,*p_init_var,*p_idx);
	  }
      }
     
    }cpd;
    cpd.p_idx=&idx;
    cpd.p_init_var=&init_var;
    cpd.ppd=&pd;
    T_var1 xprev=0.;
    std::vector<T_var1> xsamp(sample_cnt);
    
    xprev=get_element(init_var,idx);
    arms_simple(cpd,xprev,xsamp,dometrop,urand);
    set_element(init_var,idx,xsamp.back());
  }
  
};

#endif
