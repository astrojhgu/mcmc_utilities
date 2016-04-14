#include "gibbs_sampler.hpp"
#include "distribution.hpp"
#include <vector>
#include <cassert>
#include <algorithm>

namespace mcmc_utilities
{
  /**
     ``particles'' in particle filter
     An abstraction of system state
     We use templated class to code the particle class,
     so that it can be applied a wide range of problems.
     T_p represents the numerical type of probability
     T_status represnets the state vector of the system.
   */
  template <typename T_p,typename T_state>
  struct particle
  {
    /**
       Member variables:
       state: system state vector
       weight: weight of each particle in the particle filter
     */
    T_state state;
    T_p weight;

    /**
       Default constructor
     */
    particle()
    {}

    /**
       Initialize the particle with a certain state vector and its weight.
     */
    particle(const T_state& s,const T_p& w)
      :state(s),weight(w)
    {}

    /**
       Copy constructor
     */
    particle(const particle& rhs)
      :state(rhs.state),weight(rhs.weight)
    {}


    /**
       Assignment operator
     */
    particle& operator=(const particle& rhs)
    {
      resize(state,get_size(rhs.state));
      state=rhs.state;
      return *this;
    }
  };
  

  /**
     system model for the particle filter
     template parameters:
     T_p: probability data type
     T_state: system state vector type
     T_obs: observed quantity data type
     T_t: system evolution index data type, which is usually time
   */
  template <typename T_p,typename T_state,typename T_obs,typename T_t>
  class pf_model
  {
  public:
    typedef typename element_type_trait<T_state>::element_type T_state1;
  private:
    /**
       No copy constructor
     */
    pf_model(const pf_model&);

    /**
       A parameter affects the important sampling
     */
    T_p alpha;
    
  public:
    typedef typename element_type_trait<T_state>::element_type stat_element_type;
    typedef typename element_type_trait<T_obs>::element_type obs_element_type;
    
    
  public:
    /**
       default constructor
     */
    pf_model()
      :alpha(.5)
    {}

    /**
       destructor
     */
    virtual ~pf_model(){}

    
    /**
       set the value of alpha
     */
    void set_proposed_distribution_factor(const T_p& p)
    {
      alpha=p;
    }

    /**
       returns log(p(x_n|x_{n-1}))
     */
    T_p evol_log_prob(const T_state& x,const T_t& t,const T_state& prev_state,const T_t& prev_t,int n)const
    {
      return do_evol_log_prob(x,t,prev_state,prev_t,n);
    }

    /**
       returns log(p(obsx|x_n))
     */
    T_p obs_log_prob(const T_obs& y,const T_state& x,const T_t& t,int n)const
    {
      return do_obs_log_prob(y,x,t,n);
    }

    /**
       returns posterior probability
     */
    T_p combined_log_prob(const T_obs& y,const T_state& x,const T_t& t,const T_state& prev_state,const T_t& prev_t,int n)const
    {
      T_p result=evol_log_prob(x,t,prev_state,prev_t,n)+obs_log_prob(y,x,t,n);
      if(!std::isfinite(result))
	{
	  std::cerr<<"inside pf-combined_log_prob"<<std::endl;
	  std::cerr<<"n="<<n<<std::endl;
	  for(int i=0;i<x.size();++i)
	    {
	      std::cerr<<"x["<<i<<"]="<<x[i]<<std::endl;
	    }
	}

      return result;
    }

    /**
       proposed important sampling distribution, which is constructed as
       (p(x_n|x_{n-1})*p(obsx|x_n))^alpha
     */
    T_p proposed_log_prob(const T_obs& y,const T_state& x,const T_t& t,const T_state& prev_state,const T_t& prev_t,int n)const
    {
      T_p result= alpha*combined_log_prob(y,x,t,prev_state,prev_t,n);
      if(!std::isfinite(result))
	{
	  std::cerr<<"inside pf"<<std::endl;
	  std::cerr<<"n="<<n<<std::endl;
	  for(int i=0;i<x.size();++i)
	    {
	      std::cerr<<"x["<<i<<"]="<<x[i]<<std::endl;
	    }
	}
      return result;
    }

    /**
       the support of system state variable range, which is used by the sampling algorithm
     */
    std::pair<T_state1,T_state1> state_var_range(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
    {
      return do_state_var_range(t,prev_state,prev_t,ndim);
    }

    /**
       initial points used in sampling (especially for the adaptive rejection metropolis sampling algorithm, see arms.hpp)
     */
    std::vector<T_state1> init_points(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
    {
      return do_init_points(t,prev_state,prev_t,ndim);
    }
      

  public:
    /**
       Update the system state with the sequantial important sampling algorithm,
       given the new observed quantity
     */
    void update_sir(const T_obs& y,const T_t& t,std::vector<particle<T_p,T_state> >& particle_list,T_t& prev_t, base_urand<T_p>& rnd)const
    {
      /**
	 internal class representing the distribution of system state
       */
      class cprob
	:public probability_density_md<T_p,T_state>
      {
      public:
	typedef typename element_type_trait<T_state>::element_type T_var1;
      private:
	const pf_model<T_p,T_state,T_obs,T_t>*  ptr_pf_model;
	/*
	  a pointer to the observed quantity
	**/
	const T_obs* ptr_obs_vec;
	/**
	   a pointer to the particle
	*/
	const T_state* ptr_particle;
	/**
	   a pointer to the previous time value
	 */
	const T_t* ptr_prev_t;
	/**
	   a pointer to the current time value
	 */
	const T_t* ptr_t;
	
	friend class pf_model<T_p,T_state,T_obs,T_t>;
      public:
	/**
	   joint distribution of the system state
	 */
	T_p do_eval_log(const T_state& x,int n)const
	{
	  //return ptr_pf_model->evol_log_prob(x,*ptr_t,*ptr_particle,*ptr_prev_t);
	  //return ptr_pf_model->evol_log_prob(x,*ptr_t,*ptr_particle,*ptr_prev_t);
	  T_p result= ptr_pf_model->proposed_log_prob(*ptr_obs_vec,x,*ptr_t,*ptr_particle,*ptr_prev_t,n);
	  if(!std::isfinite(result))
	    {
	      std::cerr<<"inside pf-prob:"<<std::endl;
	      std::cerr<<"n="<<n<<std::endl;
	      for(int i=0;i<x.size();++i)
		{
		  std::cerr<<"x["<<i<<"]="<<x[i]<<std::endl;
		}
	    }
	  return result;
	}

	/**
	   system state variable range
	 */
	std::pair<T_var1,T_var1> do_var_range(const T_state& x,size_t ndim)const
	{
	  //ptr_pf_model->stat_var_range(x0,x1,x2);
	  //ptr_pf_model->stat_var_range(xl,xr,*ptr_particle,ndim);
	  return ptr_pf_model->state_var_range(*ptr_t,*ptr_particle,*ptr_prev_t,ndim);
	}

	/**
	   initial points for adaptive rejection metropolis sampling method
	 */
	std::vector<T_var1> do_init_points(const T_state& x,size_t ndim)const
	{
	  return ptr_pf_model->init_points(*ptr_t,*ptr_particle,*ptr_prev_t,ndim);
	}
	
      };

      /**
	 cumulative distribution of particles according to the particle weights
       */
      std::vector<T_p> weight_cdf(particle_list.size());
      
      std::vector<particle<T_p,T_state> > updated_state(particle_list.size());
      
      std::vector<T_p> log_weight(particle_list.size());

      /**
	 to determine whether the random number generator is parallel
       */
      if(rnd.is_parallel())
	{
#pragma omp parallel for
	  for(size_t i=0;i<particle_list.size();++i)
	    {
	      cprob prob;
	      /**
		 assign the pointer members in the cprob
	       */
	      
	      prob.ptr_pf_model=this;
	      prob.ptr_obs_vec=std::addressof(y);
	      prob.ptr_particle=std::addressof(particle_list[i].state);
	      prob.ptr_prev_t=std::addressof(prev_t);
	      prob.ptr_t=std::addressof(t);
	      T_state new_pred(particle_list[i].state);
	      //ofstream ofs("log.txt");
	      /**
		 perform the gibbs sampling according to the joint distribution
	       */
	      gibbs_sample(prob,new_pred,rnd);

	      /**
		 calculate the weights
	       */
	      log_weight[i]=(1-alpha)*combined_log_prob(y,new_pred,t,particle_list[i].state,prev_t,-1);

	      /**
		 assign the predicted particles to the particle list
	       */
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
	      log_weight[i]=(1-alpha)*combined_log_prob(y,new_pred,t,particle_list[i].state,prev_t,-1);
	      particle_list[i].state=new_pred;
	    }
	}
      
      T_p max_log_weight=*(std::max_element(log_weight.begin(),log_weight.end()));
      /**
	 assign the new weight to the particles in the list
       */
#pragma omp parallel for
      for(size_t i=0;i<particle_list.size();++i)
	{
	  particle_list[i].weight=std::exp(log_weight[i]-max_log_weight);
	}

      /**
	 calculate the particle cumulated distribution
       */
      weight_cdf[0]=particle_list[0].weight;
      for(size_t i=1;i<particle_list.size();++i)
	{
	  weight_cdf[i]=weight_cdf[i-1]+particle_list[i].weight;
	}

      /**
	 perform the resampling according the accumulate distribution
       */
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
	  /**
	     another non-parallelized version when the random number generator is not parallel
	   */
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
      
      /**
	 swap the new system state with the old ones
       */
      prev_t=t;
      particle_list.swap(updated_state);
      
    }
  private:
    /**
       following virtual members should be implemented when applying to a certain problem
     */

    /**
       system evolution model
     */
    virtual T_p do_evol_log_prob(const T_state& x,const T_t& t,const T_state& prev_stat,const T_t& prev_t,int n)const=0;

    /**
       observation model
     */
    virtual T_p do_obs_log_prob(const T_obs& y,const T_state& x,const T_t& t,int n)const=0;

    /**
       system state variable range
     */
    virtual std::pair<T_state1,T_state1> do_state_var_range(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const=0;

    /**
       initial points, which has an default implement as follows.
       override it if necessary.
     */
    virtual std::vector<T_state1> do_init_points(const T_t& t,const T_state& prev_state,const T_t& prev_t,size_t ndim)const
    {
      return std::vector<T_state1>();
    }
  };
};
