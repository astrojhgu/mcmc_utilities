#ifndef GEN_GAMMA_NODE_HPP
#define GEN_GAMMA_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <cmath>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class gen_gamma_node
    :public stochastic_node<T,T_vector>
  {
  private:
  public:
    gen_gamma_node()
      :stochastic_node<T,T_vector>(3,1)
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      T r=this->parent(0);
      T lambda=this->parent(1);
      T b=this->parent(2);
      return std::log(b)+b*r*std::log(lambda)+(b*r-1)*std::log(x)-std::pow(lambda*x,b)-std::lgamma(r);
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
      auto p=new gen_gamma_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };
  
  template <typename T,template <typename TE> class T_vector>
  class gen_gamma_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    gen_gamma_node_factory()
      :abstract_node_factory<T,T_vector>({"r","lambda","b"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new gen_gamma_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
}

#endif
