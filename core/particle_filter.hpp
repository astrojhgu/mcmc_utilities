#include "gibbs_sampler.hpp"
#include "distribution.hpp"
#include <vector>
#include <cassert>
#include <algorithm>

namespace mcmc_utilities
{
  template <typename T_p,typename T_stat>
  struct particle
  {
    T_stat state;
    T_p weight;

    particle()
    {}

    particle(const particle& rhs)
      :state(rhs.state),weight(rhs.weight)
    {}

    particle& operator=(const particle& rhs)
    {
      resize(state,get_size(rhs.state));
      state=rhs.state;
      return *this;
    }
  };
  
 
  template <typename T_p,typename T_stat,typename T_obs,typename T_t,typename T_urand>
  class pf_model
  {
  public:
    typedef typename element_type_trait<T_stat>::element_type T_stat1;
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
    
    void stat_var_range(T_stat1& xl,T_stat1& xr,const T_stat& x0,size_t ndim)const
    {
      return do_stat_var_range(xl,xr,x0,ndim);
    }

  public:
    void update_sir(const T_obs& y,const T_t& t,std::vector<particle<T_p,T_stat> >& particle_list,T_t& prev_t, const T_urand& urand)const
    {
      class cprob
	:public probability_density_md<T_p,T_stat>
      {
      public:
	typedef typename element_type_trait<T_stat>::element_type T_var1;
      private:
	const pf_model<T_p,T_stat,T_obs,T_t, T_urand>*  ptr_pf_model;
	const T_obs* ptr_obs_vec;
	const T_stat* ptr_particle;
	const T_t* ptr_prev_t;
	const T_t* ptr_t;
	friend class pf_model<T_p,T_stat,T_obs,T_t,T_urand>;
      public:
	T_p do_eval_log(const T_stat& x)const
	{
	  return ptr_pf_model->evol_log_prob(x,*ptr_t,*ptr_particle,*ptr_prev_t);
	}
	void do_var_range(T_var1& xl,T_var1& xr,const T_stat& x,size_t ndim)const
	{
	  //ptr_pf_model->stat_var_range(x0,x1,x2);
	  //ptr_pf_model->stat_var_range(xl,xr,*ptr_particle,ndim);
	  ptr_pf_model->stat_var_range(xl,xr,x,ndim);
	}
      };

      std::vector<T_p> weight_cdf(particle_list.size());
      std::vector<particle<T_p,T_stat> > updated_stat(particle_list.size());
      std::vector<T_p> log_weight(particle_list.size());
#pragma omp parallel for
      for(int i=0;i<particle_list.size();++i)
	{
	  cprob prob;
	  prob.ptr_pf_model=this;
	  prob.ptr_obs_vec=&y;
	  prob.ptr_particle=&(particle_list[i].state);
	  prob.ptr_prev_t=&prev_t;
	  prob.ptr_t=&t;
	  T_stat new_pred(particle_list[i].state);
	  //ofstream ofs("log.txt");

	  gibbs_sample(prob,new_pred, 1 ,urand);
	  
	  particle_list[i].state=new_pred;
	  //particle_list[i].weight=std::exp(obs_log_prob(y,new_pred,t));
	  log_weight[i]=obs_log_prob(y,new_pred,t);
	  //cout<<new_pred[0]<<" "<<particle_list[i].weight<<endl;
	}
      T_p max_log_weight=*(std::max_element(log_weight.begin(),log_weight.end()));
#pragma omp parallel for
      for(int i=0;i<particle_list.size();++i)
	{
	  particle_list[i].weight=std::exp(log_weight[i]-max_log_weight);
	}
      weight_cdf[0]=particle_list[0].weight;
      for(int i=1;i<particle_list.size();++i)
	{
	  weight_cdf[i]=weight_cdf[i-1]+particle_list[i].weight;
	}
#pragma omp parallel for
      for(int i=0;i<particle_list.size();++i)
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
	  updated_stat[i]=particle_list[n];
	  updated_stat[i].weight=1;
	}
      prev_t=t;
      particle_list.swap(updated_stat);
      
    }
  private:
    virtual T_p do_evol_log_prob(const T_stat& x,const T_t& t,const T_stat& particle_list,const T_t& prev_t)const=0;
    virtual T_p do_obs_log_prob(const T_obs& y,const T_stat& x,const T_t& t)const=0;
    virtual void do_stat_var_range(T_stat1& xl,T_stat1& xr,const T_stat& x0,size_t ndim)const=0;
  };
};
