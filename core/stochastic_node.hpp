#ifndef STOCHASTIC_NODE
#define STOCHASTIC_NODE
#include <memory>
#include <vector>
#include <future>
#include <numeric>
#include <map>
#include "error_handler.hpp"
#include "mcmc_traits.hpp"
#include "arms.hpp"
#include "slicer.hpp"
#include "base_urand.hpp"
#include "discrete_sample.hpp"
#include "continuous_sample.hpp"
#include "node.hpp"
#include "ptransform.hpp"
namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class stochastic_node
    :public node<T,T_vector>
  {
  private:
    T_vector<T> v;//store current value(s) of this node
    T_vector<int> observed;
    size_t current_idx;
  public:
    stochastic_node(size_t nparents,const T_vector<T>& v_)
      :node<T,T_vector>(nparents,get_size(v_)),
      v{v_},observed(get_size(v_)),current_idx(0)
    {
    }

    stochastic_node(size_t nparents,T v_)
      :node<T,T_vector>(nparents,1),
      v(1),observed(get_size(v)),current_idx(0)
    {
      //v[0]=v_;
      set_element(v,0,v_);
    }
    

    stochastic_node()=delete;
    stochastic_node(const stochastic_node<T,T_vector>& )=delete;
    stochastic_node<T,T_vector>& operator=(const stochastic_node<T,T_vector>&)=delete;
  public:
    bool is_observed(size_t n)const
    {
      //return observed[n]!=0;
      return get_element(observed,n)!=0;
    }

    void set_observed(size_t n,bool b)
    {
      //observed[n]=b?1:0;
      set_element(observed,n,b?1:0);
    }

    size_t num_of_observed()const
    {
      return std::accumulate(std::begin(observed),
			     std::end(observed),
			     0,
			     [](int x1,int x2){return x1+x2;});
    }

    size_t num_of_unobserved()const
    {
      return this->num_of_dims()-num_of_observed();
    }
    
    T log_posterior_prob()const
    {
      return log_prob()+this->log_likelihood();
    }

    T eval_log(const T& x)const
    {
      set_element(const_cast<stochastic_node*>(this)->v,current_idx,x);
      return log_posterior_prob();
    }

    void set_current_idx(size_t i)
    {
      current_idx=i;
    }

    size_t get_current_idx()const
    {
      return current_idx;
    }

    T log_likelihood()const
    {
      T result=static_cast<T>(0);
      auto& ss=this->get_all_stochastic_children();
      T_vector<T> results(get_size(ss));
      ////followings are the current implementation of computing the log likelihood
      //
#if defined USE_OMP_TRANSFORM
      pvec_transform_omp(ss,results,[](const stochastic_node<T,T_vector>* const& p){return p->log_prob();},10);      
#elif defined USE_THREAD_TRANSFORM
      pvec_transform_thread(ss,results,[](const stochastic_node<T,T_vector>* const& p){return p->log_prob();},10);
#else
      std::transform(std::begin(ss),std::end(ss),std::begin(results),[](const stochastic_node<T,T_vector>* const& p){return p->log_prob();});
#endif
      result=std::accumulate(std::begin(results),std::end(results),static_cast<T>(0));
      return result;
    }
    
    T log_prob()const
    {
      return do_log_prob();
    }

    void set_value(size_t idx,const T& v_)
    {
      this->set_initialized(idx,true);
      set_element(v,idx,this->regulate(idx,v_));
    }

    void set_observed_value(size_t idx,const T& v_)
    {
      this->set_value(idx,v_);
      this->set_observed(idx,true);
    }

    T regulate(size_t idx,const T& x)const
    {
      return do_regulate(idx,x);
    }

    std::pair<T,T> var_range()const
    {
      return this->do_var_range();
    }

  public:
    void sample(base_urand<T>& urand)
    {
      do_sample(urand);
    }

  public:
    void default_sample(base_urand<T>& urand)
    {
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  if(is_observed(i))
	    {
	      continue;
	    }
	  this->set_current_idx(i);
	  T xprev=this->value(i);
	  std::pair<T,T> xrange(this->do_var_range());
	  
	  //xprev=xprev<xrange.first?xrange.first:xprev;
	  //xprev=xprev>xrange.second?xrange.second:xprev;
	  
	  if(xprev<xrange.first||xprev>xrange.second)
	    {
	      this->init_value(i);
	      xprev=this->value(i);
	    }
	  
	  //arms_simple(*this,xprev,xsamp,dometrop(),rnd);
	  
	  if(is_continuous(i))
	    {
	      T_vector<T> init_x(this->do_init_points());
	      xprev=continuous_sample([&](const T& x){return this->eval_log(x);},xrange,init_x, xprev,1,urand);
	    }
	  else
	    {
	      T_vector<T> candidate_points(this->do_candidate_points());
	      xprev=discrete_sample([&](const T& x){return this->eval_log(x);},xrange, candidate_points, xprev,10,urand);
	    }
	  
	  this->set_value(i,xprev);
	}
    }

  public:
    virtual bool is_continuous(size_t idx)const=0;
    T_vector<T> init_points()const
    {
      return do_init_points();
    }

    T_vector<T> candidate_points()const
    {
      return this->do_candidate_points();
    }
    
  private:
    virtual void do_sample(base_urand<T>& urand)
    {
      default_sample(urand);
    }
    
    virtual T do_log_prob()const=0;
    T do_value(size_t idx)const override final
    {
      //return v[idx];
      return get_element(v,idx);
    }

    
    void do_connect_to_parent(node<T,T_vector>*  rhs,size_t n,size_t idx) override
    {
      //set_element(this->parents,n,std::make_pair(rhs,idx));
      this->set_parent(n,std::make_pair(rhs,idx));
      rhs->add_stochastic_child(this);
    }

    virtual T do_regulate(size_t idx,const T& x)const
    {
      return x;
    }

    virtual std::pair<T,T> do_var_range()const=0;

    virtual T_vector<T> do_init_points()const
    {
      T_vector<T> result(5);
      std::pair<T,T> xrange(do_var_range());
      T xl=xrange.first,xr=xrange.second;
      for(size_t n=0;n<get_size(result);++n)
	{	 
	  set_element(result,n,xl+(xr-xl)/(get_size(result)+1)*(n+1));
	}
      return result;
    }

    virtual T_vector<T> do_candidate_points()const
    {
      return T_vector<T>();
    }
  };  
}


#endif
