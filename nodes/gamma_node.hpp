#ifndef GAMMA_NODE_HPP
#define GAMMA_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <cmath>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class gamma_node
    :public stochastic_node<T,T_vector>
  {
  private:
  public:
    gamma_node()
      :stochastic_node<T,T_vector>(2,0)
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      T r=this->parent(0);
      T lambda=this->parent(1);
      return r*std::log(lambda)+(r-1)*std::log(x)-lambda*x-std::lgamma(r);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      return std::make_pair((T)1e-10/this->parent(1),(T)(20.)/this->parent(1));
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,this->parent(0)/this->parent(1));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new gamma_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };
  
  template <typename T,template <typename TE> class T_vector>
  class gamma_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    gamma_node_factory()
      :abstract_node_factory<T,T_vector>({"r","lambda"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new gamma_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
}

#endif
