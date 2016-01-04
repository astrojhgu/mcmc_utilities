#include "gibbs_sampler.hpp"
#include "distribution.hpp"
#include <vector>
#include <cassert>
#include <algorithm>

namespace mcmc_utilities
{
  template <typename T_p,typename T_state>
  struct particle
  {
    T_state state;
    T_p weight;

    particle()
    {}

    particle(const T_state& s,const T_p& w)
      :state(s),weight(w)
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
  
 
  template <typename T_p,typename T_state,typename T_obs,typename T_t>
  class pf_model
  {
  public:
    typedef typename element_type_trait<T_state>::element_type T_state1;
  private:
    pf_model(const pf_model&);
    T_p alpha;
    
  public:
    typedef typename element_type_trait<T_state>::element_type stat_element_type;
    typedef typename element_type_trait<T_obs>::element_type obs_element_type;
    
    
  public:
    pf_model()
      :alpha(.5)
    {}
    virtual ~pf_model(){}

    void set_proposed_distribution_factor(const T_p& p)
    {
      alpha=p;
    }
    
    T_p evol_log_prob(const T_state& x,const T_t& t,const T_state& prev_state,const T_t& prev_t)const
    {
      return do_evol_log_prob(x,t,prev_state,prev_t);
    }

    T_p obs_log_prob(const T_obs& y,const T_state& x,const T_t& t)const
    {
      return do_obs_log_prob(y,x,t);
    }

    T_p combined_log_prob(const T_obs& y,const T_state& x,const T_t& t,const T_state& prev_state,const T_t& prev_t)const
    {
      return evol_log_prob(x,t,prev_state,prev_t)+obs_log_prob(y,x,t);
    }

    T_p proposed_log_prob(const T_obs& y,const T_state& x,const T_t& t,const T_state& prev_state,const T_t& prev_t)const
    {
      return alpha*combined_log_prob(y,x,t,prev_state,prev_t);
    }
    
    std::pair<T_state1,T_state1> state_var_range(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
    {
      return do_state_var_range(t,prev_state,prev_t,ndim);
    }

    std::vector<T_state1> init_points(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
    {
      return do_init_points(t,prev_state,prev_t,ndim);
    }
      

  public:
    void update_sir(const T_obs& y,const T_t& t,std::vector<particle<T_p,T_state> >& particle_list,T_t& prev_t, base_urand<T_p>& rnd)const
    {
      class cprob
	:public probability_density_md<T_p,T_state>
      {
      public:
	typedef typename element_type_trait<T_state>::element_type T_var1;
      private:
	const pf_model<T_p,T_state,T_obs,T_t>*  ptr_pf_model;
	const T_obs* ptr_obs_vec;
	const T_state* ptr_particle;
	const T_t* ptr_prev_t;
	const T_t* ptr_t;
	friend class pf_model<T_p,T_state,T_obs,T_t>;
      public:
	T_p do_eval_log(const T_state& x)const
	{
	  //return ptr_pf_model->evol_log_prob(x,*ptr_t,*ptr_particle,*ptr_prev_t);
	  //return ptr_pf_model->evol_log_prob(x,*ptr_t,*ptr_particle,*ptr_prev_t);
	  return ptr_pf_model->proposed_log_prob(*ptr_obs_vec,x,*ptr_t,*ptr_particle,*ptr_prev_t);
	}
	std::pair<T_var1,T_var1> do_var_range(const T_state& x,size_t ndim)const
	{
	  //ptr_pf_model->stat_var_range(x0,x1,x2);
	  //ptr_pf_model->stat_var_range(xl,xr,*ptr_particle,ndim);
	  return ptr_pf_model->state_var_range(*ptr_t,*ptr_particle,*ptr_prev_t,ndim);
	}

	std::vector<T_var1> do_init_points(const T_state& x,size_t ndim)const
	{
	  return ptr_pf_model->init_points(*ptr_t,*ptr_particle,*ptr_prev_t,ndim);
	}
	
      };

      std::vector<T_p> weight_cdf(particle_list.size());
      std::vector<particle<T_p,T_state> > updated_state(particle_list.size());
      std::vector<T_p> log_weight(particle_list.size());

      if(rnd.is_parallel())
	{
#pragma omp parallel for
	  for(size_t i=0;i<particle_list.size();++i)
	    {
	      cprob prob;
	      prob.ptr_pf_model=this;
	      prob.ptr_obs_vec=std::addressof(y);
	      prob.ptr_particle=std::addressof(particle_list[i].state);
	      prob.ptr_prev_t=std::addressof(prev_t);
	      prob.ptr_t=std::addressof(t);
	      T_state new_pred(particle_list[i].state);
	      //ofstream ofs("log.txt");
	      
	      gibbs_sample(prob,new_pred,rnd);
	      log_weight[i]=(1-alpha)*combined_log_prob(y,new_pred,t,particle_list[i].state,prev_t);
	      particle_list[i].state=new_pred;
	    }
	}
      else
	{
	  for(size_t i=0;i<particle_list.size();++i)
	    {
	      cprob prob;
	      prob.ptr_pf_model=this;
	      prob.ptr_obs_vec=std::addressof(y);
	      prob.ptr_particle=std::addressof(particle_list[i].state);
	      prob.ptr_prev_t=std::addressof(prev_t);
	      prob.ptr_t=std::addressof(t);
	      T_state new_pred(particle_list[i].state);
	      //ofstream ofs("log.txt");
	      
	      gibbs_sample(prob,new_pred,rnd);
	      log_weight[i]=(1-alpha)*combined_log_prob(y,new_pred,t,particle_list[i].state,prev_t);
	      particle_list[i].state=new_pred;
	    }
	}
      
      T_p max_log_weight=*(std::max_element(log_weight.begin(),log_weight.end()));
#pragma omp parallel for
      for(size_t i=0;i<particle_list.size();++i)
	{
	  particle_list[i].weight=std::exp(log_weight[i]-max_log_weight);
	}
      weight_cdf[0]=particle_list[0].weight;
      for(size_t i=1;i<particle_list.size();++i)
	{
	  weight_cdf[i]=weight_cdf[i-1]+particle_list[i].weight;
	}

      if(rnd.is_parallel())
	{
#pragma omp parallel for
	  for(size_t i=0;i<particle_list.size();++i)
	    {
	      T_p p=0;
	      p=rnd()*weight_cdf.back();
	      size_t n=0;
	      bool changed=false;
	      for(size_t j=0;j<weight_cdf.size();++j)
		{
		  if(weight_cdf[j]>=p)
		    {
		      n=j;
		      changed=true;
		      break;
		    }
		}
	      assert(changed);
	      updated_state[i]=particle_list[n];
	      updated_state[i].weight=1;
	    }
	}
      else
	{
	  for(size_t i=0;i<particle_list.size();++i)
	    {
	      T_p p=0;
	      p=rnd()*weight_cdf.back();
	      size_t n=0;
	      bool changed=false;
	      for(size_t j=0;j<weight_cdf.size();++j)
		{
		  if(weight_cdf[j]>=p)
		    {
		      n=j;
		      changed=true;
		      break;
		    }
		}
	      assert(changed);
	      updated_state[i]=particle_list[n];
	      updated_state[i].weight=1;
	    }	  
	}
      
      
      prev_t=t;
      particle_list.swap(updated_state);
      
    }
  private:
    virtual T_p do_evol_log_prob(const T_state& x,const T_t& t,const T_state& prev_stat,const T_t& prev_t)const=0;
    virtual T_p do_obs_log_prob(const T_obs& y,const T_state& x,const T_t& t)const=0;
    virtual std::pair<T_state1,T_state1> do_state_var_range(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const=0;
    
    virtual std::vector<T_state1> do_init_points(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
    {
      return std::vector<T_state1>();
    }
  };
};
