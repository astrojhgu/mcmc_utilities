#include "gibbs_sampler.hpp"
#include "distribution.hpp"
#include <vector>
#include <cassert>


namespace mcmc_utilities
{
  template <typename T_p,typename T_stat>
  struct particle
  {
    T_stat status;
    T_p weight;
  };
  
 
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
    T_p evol_log_prob(const T_stat& x,const T_t& t,const T_stat& prev_stat,const T_t& prev_t)const
    {
      return do_evol_log_prob(x,t,prev_stat,prev_t);
    }

    T_p obs_log_prob(const T_obs& y,const T_stat& x,const T_t& t)const
    {
      return do_obs_log_prob(y,x,t);
    }

    T_p combined_log_prob(const T_obs& y,const T_stat& x,const T_t& t,const T_stat& prev_stat,const T_t& prev_t)const
    {
      return evol_log_prob(x,t,prev_stat,prev_t)+obs_log_prob(y,x,t);
    }
    
    void stat_var_range(const T_stat& x0,T_stat& xl,T_stat& xr)const
    {
      return do_stat_var_range(x0,xl,xr);
    }

  public:
    void update_sir(const T_obs& y,const T_t& t,std::vector<particle<T_p,T_stat> >& prev_stat,T_t& prev_t)const
    {
      class cprob
	:public probability_density_md<T_p,T_stat>
      {
      private:
	const pf_model*  ptr_pf_model;
	const T_obs* ptr_obs_vec;
	const T_stat* ptr_prev_stat;
	const T_t* ptr_prev_t;
	const T_t* ptr_t;
	friend class pf_model;
      public:
	T_p do_eval_log(const T_stat& x)const
	{
	  //return ptr_pf_model->evol_log_prob(x,*ptr_t,*ptr_prev_stat,*ptr_prev_t)
	  //+ptr_pf_model->obs_log_prob(*ptr_obs_vec,x,*ptr_t);
	  //return ptr_pf_model->combined_log_prob(*ptr_obs_vec,x,*ptr_t,*ptr_prev_stat,*ptr_prev_t);
	  return ptr_pf_model->evol_log_prob(x,*ptr_t,*ptr_prev_stat,*ptr_prev_t);
	}
	cprob* do_clone()const
	{
	  return new cprob(*this);
	}
	
	void do_var_range(T_stat& xl,T_stat& xr)const
	{
	  //ptr_pf_model->stat_var_range(x0,x1,x2);
	  ptr_pf_model->stat_var_range(*ptr_prev_stat,xl,xr);
	}
      }prob;

      std::vector<T_p> weight_cdf;
      std::vector<particle<T_p,T_stat> > updated_stat;
      for(int i=0;i<prev_stat.size();++i)
	{
	  prob.ptr_pf_model=this;
	  prob.ptr_obs_vec=&y;
	  prob.ptr_prev_stat=&(prev_stat[i].status);
	  prob.ptr_prev_t=&prev_t;
	  prob.ptr_t=&t;
	  T_stat new_pred(prev_stat[i].status);
	  //ofstream ofs("log.txt");

	  gibbs_sample(prob,new_pred);
	  
	  prev_stat[i].status=new_pred;
	  prev_t=t;
	  prev_stat[i].weight=std::exp(obs_log_prob(y,new_pred,t));
	  //cout<<new_pred[0]<<" "<<prev_stat[i].weight<<endl;
	  if(i==0)
	    {
	      weight_cdf.push_back(prev_stat[i].weight);
	    }
	  else
	    {
	      weight_cdf.push_back(weight_cdf.back()+prev_stat[i].weight);
	    }
	  //cout<<i<<" "<<weight_cdf.back()<<endl;
	}
      for(int i=0;i<prev_stat.size();++i)
	{
	  T_p p=random()/(T_p)RAND_MAX*weight_cdf.back();

	  int n=-1;
	  for(int j=0;j<weight_cdf.size();++j)
	    {
	      if(weight_cdf[j]>=p)
		{
		  n=j;
		  break;
		}
	    }
	  
	  assert(n!=-1);
	  updated_stat.push_back(prev_stat.at(n));
	  updated_stat.back().weight=1;
	}
      prev_stat.swap(updated_stat);
      
    }
  private:
    virtual T_p do_evol_log_prob(const T_stat& x,const T_t& t,const T_stat& prev_stat,const T_t& prev_t)const=0;
    virtual T_p do_obs_log_prob(const T_obs& y,const T_stat& x,const T_t& t)const=0;
    virtual void do_stat_var_range(const T_stat& x0,T_stat& xl,T_stat& xr)const=0;
  };
};
