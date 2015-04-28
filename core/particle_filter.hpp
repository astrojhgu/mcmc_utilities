#include "gibbs_sampler.hpp"
#include "distribution.hpp"
#include <vector>
#include <deque>

namespace mcmc_utilities
{
 
  template <typename T_p,typename T_stat,typename T_obs,typename T_t>
  class pf_model
  {
  private:
    pf_model(const pf_model&);
  public:
    typedef typename element_type_trait<T_stat>::element_type stat_element_type;
    typedef typename element_type_trait<T_obs>::element_type obs_element_type;
  public:
    pf_model(){}
    virtual ~pf_model(){}
    T_p evol_log_prob(const T_stat& x,const T_t& t,const std::deque<T_stat>& prev_stat,const std::deque<T_t>& prev_t)const
    {
      return do_evol_log_prob(x,t,prev_stat,prev_t);
	}

    T_p obs_log_prob(const T_obs& y,const T_stat& x,const T_t& t)const
    {
      return do_obs_log_prob(y,x,t);
    }

    void stat_var_range(T_stat& x1,T_stat& x2)const
    {
      return do_stat_var_range(x1,x2);
    }

  public:
    void update(const T_obs& y,const T_t& t,std::deque<T_stat>& prev_stat,std::deque<T_t>& prev_t)const
    {
      class cprob
	:public probability_density_md<T_p,T_stat>
      {
      private:
	const pf_model*  ptr_pf_model;
	const T_obs* ptr_obs_vec;
	const std::deque<T_stat>* ptr_prev_stat;
	const std::deque<T_t>* ptr_prev_t;
	const T_t* ptr_t;
	friend class pf_model;
      public:
	T_p do_eval_log(const T_stat& x)const
	{
	  return ptr_pf_model->evol_log_prob(x,*ptr_t,*ptr_prev_stat,*ptr_prev_t)+
	    ptr_pf_model->obs_log_prob(*ptr_obs_vec,x,*ptr_t);
	}
	cprob* do_clone()const
	{
	  return new cprob(*this);
	}
	
	void do_var_range(T_stat& x1,T_stat& x2)const
	{
	  ptr_pf_model->stat_var_range(x1,x2);
	}
      }prob;
      prob.ptr_pf_model=this;
      prob.ptr_obs_vec=&y;
      prob.ptr_prev_stat=&prev_stat;
      prob.ptr_prev_t=&prev_t;
      prob.ptr_t=&t;
      T_stat new_pred(prev_stat.back());
      gibbs_sample(prob,new_pred);
      prev_stat.push_back(new_pred);
      prev_t.push_back(t);
    }
    
  private:
    virtual T_p do_evol_log_prob(const T_stat& x,const T_t& t,const std::deque<T_stat>& prev_stat,const std::deque<T_t>& prev_t)const=0;
    virtual T_p do_obs_log_prob(const T_obs& y,const T_stat& x,const T_t& t)const=0;
    virtual void do_stat_var_range(T_stat& x1,T_stat& x2)const=0;
  };
};
