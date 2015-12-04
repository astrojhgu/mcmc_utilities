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
#include "node.hpp"

namespace mcmc_utilities
{
  template <typename T>
  class stochastic_node
    :public node<T>,public probability_density_1d<T>
  {
  private:
    std::vector<T> v;//store current value(s) of this node
    std::vector<int> observed;
    size_t current_idx;
  public:
    stochastic_node(size_t nparents,const std::vector<T>& v_)
      :node<T>(nparents,v_.size()),
      v{v_},observed(v_.size()),current_idx(0)
    {}

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
    bool is_observed(size_t n)
    {
      return observed[n]!=0;
    }

    void set_observed(size_t n,bool b)
    {
      observed[n]=b?1:0;
    }
    
    T log_posterior_prob()const
    {
      return log_prob()+this->log_likelihood();
    }

    T do_eval_log(const T& x)const override
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
      v[idx]=this->regulate(v_);
    }

  public:
    void sample(const base_urand<T>& urand)
    {
      do_sample(urand);
    }

  private:
    virtual void do_sample(const base_urand<T>& urand)
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
	  std::pair<T,T> xrange(this->var_range());
	  //xprev=xprev<xrange.first?xrange.first:xprev;
	  //xprev=xprev>xrange.second?xrange.second:xprev;

	  if(xprev<xrange.first||xprev>xrange.second)
	    {
	      this->initialize(i);
	    }
	  
	  //arms_simple(*this,xprev,xsamp,dometrop(),rnd);
	  
	  if(is_continuous(i))
	    {
	      xprev=arms(*this,xprev,10,urand);
	    }
	  else
	    {
	      xprev=discrete_sample(*this,xprev,urand);
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

    virtual bool is_continuous(size_t idx)const=0;
  };  
}


#endif
