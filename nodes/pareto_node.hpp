#ifndef PARETO_NODE_HPP
#define PARETO_NODE_HPP
#include <core/forward_sampleable_node.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>
#include <math/distributions.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class pareto_node
    :public forward_sampleable_node<T,T_vector>
  {
  private:
  public:
    pareto_node()
      :forward_sampleable_node<T,T_vector>(2,0)
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      T c=this->parent(0);
      T alpha=this->parent(1);
      return logdpar<T,T>(x,c,alpha);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      T c=this->parent(0);
      T alpha=this->parent(1);
      T eta=1e-4;
      return std::make_pair(c,c*std::pow(eta,-(T)(1.)/alpha));
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,this->parent(0));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new pareto_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };
  
  template <typename T,template <typename TE> class T_vector>
  class pareto_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    pareto_node_factory()
      :abstract_node_factory<T,T_vector>({"mu","sigma"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new pareto_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
}

#endif
