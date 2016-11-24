#ifndef FORWARD_SAMPLEABLE_NODE
#define FORWARD_SAMPLEABLE_NODE
#include "stochastic_node.hpp"
#include <algorithm>
namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class forward_sampleable_node
    :public stochastic_node<T,T_vector>
  {
  private:
    bool no_observed_children;
  public:
    forward_sampleable_node(size_t nparents,const T_vector<T>& v_)
      :stochastic_node<T,T_vector>(nparents,v_),no_observed_children(true)
    {}

    forward_sampleable_node(size_t nparents,T v_)
      :stochastic_node<T,T_vector>(nparents,v_),no_observed_children(true)
    {}

  public:
    void do_sample(base_urand<T>& urand)override
    {
      if(no_observed_children)
	{
	  //std::cerr<<"a"<<std::endl;
	  forward_sample(urand);
	}
      else
	{
	  this->default_sample(urand);
	}
    }

    virtual void forward_sample(base_urand<T>& urand)
    {
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  if(this->is_observed(i))
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
	      this->init_value(i);
	      xprev=this->value(i);
	    }
	  
	  //arms_simple(*this,xprev,xsamp,dometrop(),rnd);
	  
	  if(this->is_continuous(i))
	    {
	      T_vector<T> init_x(this->init_points());
	      xprev=continuous_sample([&](const T& x){return this->eval_log_prior(x);},xrange,init_x, xprev,1,urand);
	    }
	  else
	    {
	      T_vector<T> candidate_points(this->candidate_points());
	      xprev=discrete_sample([&](const T& x){return this->eval_log_prior(x);},xrange, candidate_points, xprev,10,urand);
	    }
	  
	  this->set_value(i,xprev);
	}
    }

    T eval_log_prior(const T& x)const
    {
      //set_element(const_cast<forward_sampleable_node<T,T_vector>*>(this)->v,this->get_current_idx(),x);
      const_cast<forward_sampleable_node<T,T_vector>*>(this)->set_value(this->get_current_idx(),x);
      return this->log_prob();
    }

    bool has_observed_children(const node<T,T_vector>* pn)const
    {
      bool result=false;
      for(const auto& p:pn->get_stochastic_children())
	{
	  result=result||has_observed_children(p);
	  if(result)
	    {
	      return result;
	    }
	  for(size_t i=0;i<p->num_of_dims();++i)
	    {
	      result=result||p->is_observed(i);
	      if(result)
		{
		  return result;
		}
	    }
	}
      for(const auto& p:pn->get_deterministic_children())
	{
	  result=result||has_observed_children(p);
	  if(result)
	    {
	      return result;
	    }
	}
      return result;
    }

    void do_freeze_topology()override
    {
      no_observed_children=!has_observed_children(this);
    }
  };
}

#endif
