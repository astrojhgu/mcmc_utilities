#ifndef STOCHASTIC_NODE
#define STOCHASTIC_NODE
#include <memory>
#include <vector>
#include <map>
#include "error_handler.hpp"
#include "arms.hpp"
#include "slicer.hpp"
#include "base_urand.hpp"
#include "discrete_sample.hpp"
#include "continuous_sample.hpp"
#include "node.hpp"

namespace mcmc_utilities
{
  template <typename T>
  class stochastic_node
    :public node<T>
  {
  private:
    std::vector<T> v;//store current value(s) of this node
    std::vector<int> observed;
    size_t current_idx;
  public:
    stochastic_node(size_t nparents,const std::vector<T>& v_)
      :node<T>(nparents,v_.size()),
      v{v_},observed(v_.size()),current_idx(0)
    {
    }

    stochastic_node(size_t nparents,T v_)
      :node<T>(nparents,1),
      v(1),observed(v.size()),current_idx(0)
    {
      v[0]=v_;
    }
    

    stochastic_node()=delete;
    stochastic_node(const stochastic_node<T>& )=delete;
    stochastic_node<T>& operator=(const stochastic_node<T>&)=delete;
  public:
    bool is_observed(size_t n)const
    {
      return observed[n]!=0;
    }

    void set_observed(size_t n,bool b)
    {
      observed[n]=b?1:0;
    }

    size_t num_of_observed()const
    {
      return std::accumulate(observed.begin(),
			     observed.end(),
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
      const_cast<stochastic_node*>(this)->v[current_idx]=x;
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
    
    T log_prob()const
    {
      return do_log_prob();
    }

    void set_value(size_t idx,const T& v_)
    {
      this->set_initialized(idx,true);
      v[idx]=this->regulate(idx,v_);
    }

    T regulate(size_t idx,const T& x)const
    {
      return do_regulate(idx,x);
    }



  public:
    void sample(base_urand<T>& urand)
    {
      do_sample(urand);
    }

  private:
    virtual void do_sample(base_urand<T>& urand)
    {
      //constexpr size_t nsamp=10;
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
	      this->initialize(i);
	      xprev=this->value(i);
	    }
	  
	  //arms_simple(*this,xprev,xsamp,dometrop(),rnd);
	  
	  if(is_continuous(i))
	    {
	      std::vector<T> init_x(this->do_init_points());
	      xprev=continuous_sample([&](const T& x){return this->eval_log(x);},xrange,init_x, xprev,1,urand);
	    }
	  else
	    {
	      std::vector<T> candidate_points(this->do_candidate_points());
	      xprev=discrete_sample([&](const T& x){return this->eval_log(x);},xrange, candidate_points, xprev,10,urand);
	    }
	  
	  this->set_value(i,xprev);
	}
    }
    
    virtual T do_log_prob()const=0;
    T do_value(size_t idx)const override final
    {
      return v[idx];
    }

    
    void do_connect_to_parent(node<T>*  rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
      rhs->add_stochastic_child(this);
    }

    virtual T do_regulate(size_t idx,const T& x)const
    {
      return x;
    }

    virtual bool is_continuous(size_t idx)const=0;

    virtual std::pair<T,T> do_var_range()const=0;

    virtual std::vector<T> do_init_points()const
    {
      std::vector<T> result(5);
      std::pair<T,T> xrange(do_var_range());
      T xl=xrange.first,xr=xrange.second;
      for(size_t n=0;n<result.size();++n)
	{	 
	  result[n]=xl+(xr-xl)/(result.size()+1)*(n+1);
	}
      return result;
    }

    virtual std::vector<T> do_candidate_points()const
    {
      return std::vector<T>();
    }
  };  
}


#endif
