#ifndef GIBBS_SAMPLER_HPP
#define GIBBS_SAMPLER_HPP
//#include "rejection_sampler_1d.hpp"
#include "arms.hpp"
#include <vector>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  void gibbs_sample(const probability_density_md<T_p,T_var>& pd,T_var& init_var,int dometrop=1)
  {
    size_t idx=0;
    typedef typename element_type_trait<T_var>::element_type T_var1;
    class conditional_probability_density
      :public probability_density_md<T_p,typename element_type_trait<T_var>::element_type>
    {
    public:
      size_t* p_idx;
      T_var* p_init_var;
      const probability_density_md<T_p,T_var>* ppd;
    private:
      T_p do_eval(const T_var1& x1)const
      {
	//set_element(*p_init_var,*p_idx,x1);
	p_init_var->at(*p_idx)=x1;
	return ppd->eval(*p_init_var);
      }

      T_p do_eval_log(const T_var1& x1)const
      {
	//set_element(*p_init_var,*p_idx,x1);
	p_init_var->at(*p_idx)=x1;
	
	T_p result= ppd->eval_log(*p_init_var);
	//std::cout<<*p_idx<<" ";
	for(int i=0;i<3;++i)
	  {
	    //std::cout<<(*p_init_var)[i]<<" ";
	  }
	//std::cout<<result<<endl;
	return result;
      }
      
      probability_density_md<T_p,T_var1>* do_clone()const
      {
	return new conditional_probability_density(*this);
      }

      void do_var_range(T_var1& x1,T_var1& x2)const
      {
	T_var xmin,xmax;
	//x1=xmin[*p_idx];
	//x2=xmax[*p_idx];
	resize(xmin,get_size(*p_init_var));
	resize(xmax,get_size(*p_init_var));
	ppd->var_range(xmin,xmax);
	x1=get_element(xmin,*p_idx);
	x2=get_element(xmax,*p_idx);
      }
    }cpd;
    cpd.p_idx=&idx;
    cpd.p_init_var=&init_var;
    //T_var init_var_orig=init_var;
    cpd.ppd=&pd;
    T_var1 xprev=0.;
    for(idx=0;idx<get_size(init_var);++idx)
      {
	xprev=get_element(init_var,idx);
	std::vector<T_var1> xsamp(10);
	arms_simple(10,cpd,xprev,xsamp,dometrop);
	set_element(init_var,idx,xsamp[0]);
	//init_var.at(idx)=xsamp[0];
      }
  }
		    
};

#endif
