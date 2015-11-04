#ifndef STOCHASTIC_NODE
#define STOCHASTIC_NODE
#include <memory>
#include <vector>
#include <map>
#include "mcmc_exception.hpp"
#include "arms.hpp"
#include "base_urand.hpp"
#include "node.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class stochastic_node
    :public node<T_p,T_var1>,public probability_density_1d<T_p,T_var1>
  {
  private:
    std::vector<T_var1> v;//store current value(s) of this node
    size_t current_idx;
  public:
    stochastic_node(size_t nparents,const std::vector<T_var1>& v_)
      :node<T_p,T_var1>(nparents,v_.size()),
      v{v_},current_idx(0)
    {}

    stochastic_node(size_t nparents,T_var1 v_)
      :node<T_p,T_var1>(nparents,1),
      v(1),current_idx(0)
    {
      v[0]=v_;
    }
    

    stochastic_node()=delete;
    stochastic_node(const stochastic_node<T_p,T_var1>& )=delete;
    stochastic_node<T_p,T_var1>& operator=(const stochastic_node<T_p,T_var1>&)=delete;
  public:
    T_p log_posterior_prob()const
    {
      return log_prob()+this->log_likelihood();
    }

    T_p do_eval_log(const T_var1& x)const override
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
    
    T_p log_prob()const
    {
      return do_log_prob();
    }

    void set_value(size_t idx,const T_var1& v_)
    {
      v[idx]=v_;
    }

  public:
    void sample(const base_urand<T_p>& urand)
    {
      do_sample(urand);
    }

  private:
    virtual void do_sample(const base_urand<T_p>& urand)
    {
      //constexpr size_t nsamp=10;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  this->set_current_idx(i);
	  T_var1 xprev=this->value(i,0);
	  std::pair<T_var1,T_var1> xrange(this->var_range());
	  xprev=xprev<xrange.first?xrange.first:xprev;
	  xprev=xprev>xrange.second?xrange.second:xprev;
	  
	  //arms_simple(*this,xprev,xsamp,dometrop(),rnd);

	  if(is_continuous(i))
	    {
	      xprev=arms(*this,xprev,10,urand);
	    }
	  else
	    {
	      xprev=discrete_sample(*this,urand);
	    }   
	  
	  this->set_value(i,xprev);
	}
    }
    
    virtual T_p do_log_prob()const=0;
    T_var1 do_value(size_t idx,size_t obsid)const override final
    {
      return v[idx];
    }

    
    void do_connect_to_parent(node<T_p,T_var1>*  rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
      rhs->add_stochastic_child(this);
    }

    virtual bool is_continuous(size_t idx)const=0;
  };  
}


#endif
