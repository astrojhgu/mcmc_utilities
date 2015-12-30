#ifndef PARTICLE_SMOOTHER
#define PARTICLE_SMOOTHER

#include "particle_filter.hpp"
#include <cmath>

namespace mcmc_utilities
{
  template <typename T_p,typename T_state,typename T_obs,typename T_t>
  class particle_smoother
  {
  public:
    typedef typename element_type_trait<T_state>::element_type T_state1;
  private:
    int verbose_level;
    particle_smoother(const particle_smoother&);
  public:
    typedef typename element_type_trait<T_state>::element_type stat_element_type;
    typedef typename element_type_trait<T_obs>::element_type obs_element_type;

    
  private:
    class smoother_model
      :public pf_model<T_p,T_state,T_obs,T_t>
    {
    private:
      particle_smoother<T_p,T_state,T_obs,T_t>* pps;
      
    public:
      smoother_model(particle_smoother<T_p,T_state,T_obs,T_t>* _pps)
	:pps(_pps)
      {}

    private:
      T_p do_evol_log_prob(const T_state& x,const T_t& t,const T_state& prev_stat,const T_t& prev_t)const
      {
	return pps->evol_log_prob(x,t,prev_stat,prev_t);
      }
      
      T_p do_obs_log_prob(const T_obs& y,const T_state& x,const T_t& t)const
      {
	return pps->obs_log_prob(y,x,t);
      }

      std::pair<T_state1,T_state1> do_state_var_range(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
      {
	return pps->state_var_range(t,prev_state,prev_t,ndim);
      }
    
      std::vector<T_state1> do_init_points(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
      {
	return pps->init_points(t,prev_state,prev_t,ndim);
      }
    }sm;

    class smoother_model_reverse
      :public pf_model<T_p,T_state,T_obs,T_t>
    {
    private:
      particle_smoother<T_p,T_state,T_obs,T_t>* pps;
    public:
      smoother_model_reverse(particle_smoother<T_p,T_state,T_obs,T_t>* _pps)
	:pps(_pps)
      {}

    private:
      T_p do_evol_log_prob(const T_state& x,const T_t& t,const T_state& future_stat,const T_t& future_t)const
      {
	return pps->evol_log_prob_rev(x,t,future_stat,future_t);
      }
      
      T_p do_obs_log_prob(const T_obs& y,const T_state& x,const T_t& t)const
      {
	return pps->obs_log_prob(y,x,t);
      }

      std::pair<T_state1,T_state1> do_state_var_range(const T_t& t,const T_state& future_state,const T_t& future_t,size_t ndim)const
      {
	return pps->state_var_range(t,future_state,future_t,ndim);
      }
    
      std::vector<T_state1> do_init_points(const T_t& t,const T_state& future_state,const T_t& future_t,size_t ndim)const
      {
	return pps->init_points(t,future_state,future_t,ndim);
      }
    }sm_rev;    
  public:
    //size_t nparticles;
    std::vector<std::vector<particle<T_p,T_state> > > history;
    std::vector<T_t> t_list;
    std::vector<T_obs> obs_list;
  public:
    particle_smoother()
      :sm(this),sm_rev(this),verbose_level(0)
    {}

    virtual ~particle_smoother(){}

    void set_verbose(int n)
    {
      verbose_level=n;
    }
  public:
    void load(const std::vector<T_obs>& obs,const std::vector<T_t>& t_vec,const std::vector<particle<T_p,T_state> >& ps,T_t t0,base_urand<T_p>& rng)
    {
      history.clear();
      obs_list.clear();
      std::vector<particle<T_p,T_state> > particles(ps);
      size_t nparticles=ps.size();
      T_t prev_t(t0);
      
      for(size_t i=0;i<obs.size();++i)
	{
	  sm.update_sir(obs.at(i),t_vec.at(i),particles,prev_t,rng);
	  history.push_back(particles);
	  t_list.push_back(t_vec[i]);
	  obs_list.push_back(obs[i]);
	  if(verbose_level)
	    {
	      std::cerr<<"loaded "<<i<<" of "<<obs.size()<<std::endl;
	    }
	}
    }

    void backward_simulate(base_urand<T_p>& rng)
    {
      std::vector<particle<T_p,T_state> > particles(history.back());
      size_t nparticles=particles.size();
      T_t future_t(t_list.back());
      
      for(int i=obs_list.size()-1;i>=0;--i)
	{
	  sm_rev.update_sir(obs_list.at(i),t_list.at(i),particles,future_t,rng);
	  history[i]=particles;
	  if(verbose_level)
	    {
	      std::cerr<<"bs:"<<(obs_list.size()-1-i)<<" of "<<obs_list.size()<<std::endl;
	    }
	}
    }

    void forward_simulate(base_urand<T_p>& rng)
    {
      std::vector<particle<T_p,T_state> > particles(history.front());
      size_t nparticles=particles.size();
      T_t prev_t(t_list.front());
      
      for(int i=0;i!=obs_list.size();++i)
	{
	  sm.update_sir(obs_list.at(i),t_list.at(i),particles,prev_t,rng);
	  history[i]=particles;
	  if(verbose_level)
	    {
	      std::cerr<<"fs:"<<i<<" of "<<obs_list.size()<<std::endl;
	    }
	}
    }   
    

    std::vector<T_state> draw_realization(base_urand<T_p>& rng)
    {
      std::vector<T_state> result(history.size());
      for(int i=history.size()-1;i>=0;--i)
	{
	  std::vector<T_p> weights(history[i].size());
	  if(false&&i!=history.size()-1)
	    {
	      for(int j=0;j<history[i].size();++j)
		{
		  weights[j]=history[i][j].weight*std::exp(evol_log_prob(result[i+1],t_list[i+1],history[i][j].state,t_list[i]));
		}
	    }
	  else
	    {
	      for(int j=0;j<history[i].size();++j)
		{
		  weights[j]=history[i][j].weight;
		}
	    }
	  std::vector<T_p> weight_cdf(history[i].size());
	  weight_cdf[0]=weights[0];
	  for(size_t j=1;j<weights.size();++j)
	    {
	      weight_cdf[j]=weight_cdf[j-1]+weights[j];
	    }
	  T_p p=rng()*weight_cdf.back();
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
	  result[i]=history[i][n].state;
	}
      return result;
    }
    
  public:
    T_p evol_log_prob(const T_state& x,const T_t& t,const T_state& prev_state,const T_t& prev_t)const
    {
      return do_evol_log_prob(x,t,prev_state,prev_t);
    }

    T_p evol_log_prob_rev(const T_state& x,const T_t& t,const T_state& future_state,const T_t& future_t)const
    {
      return do_evol_log_prob_rev(x,t,future_state,future_t);
    }
    
    T_p obs_log_prob(const T_obs& y,const T_state& x,const T_t& t)const
    {
      return do_obs_log_prob(y,x,t);
    }

    std::pair<T_state1,T_state1> state_var_range(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
    {
      return do_state_var_range(t,prev_state,prev_t,ndim);
    }

    std::pair<T_state1,T_state1> state_var_range_rev(const T_t& t,const T_state& future_state,const T_t& future_t,size_t ndim)const
    {
      return do_state_var_range_rev(t,future_state,future_t,ndim);
    }    
    
    std::vector<T_state1> init_points(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
    {
      return do_init_points(t,prev_state,prev_t,ndim);
    }

  private:
    virtual T_p do_evol_log_prob(const T_state& x,const T_t& t,const T_state& prev_stat,const T_t& prev_t)const=0;
    virtual T_p do_evol_log_prob_rev(const T_state& x,const T_t& t,const T_state& prev_stat,const T_t& prev_t)const=0;
    virtual T_p do_obs_log_prob(const T_obs& y,const T_state& x,const T_t& t)const=0;
    virtual std::pair<T_state1,T_state1> do_state_var_range(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const=0;
    virtual std::pair<T_state1,T_state1> do_state_var_range_rev(const T_t& t,const T_state& future_state,const T_t& future_t,size_t ndim)const=0;
    
    virtual std::vector<T_state1> do_init_points(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
    {
      return std::vector<T_state1>();
    }
  };
}


#endif
