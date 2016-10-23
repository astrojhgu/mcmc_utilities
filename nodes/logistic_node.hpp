#ifndef LOGISTIC_NODE_HPP
#define LOGISTIC_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class logistic_node
    :public stochastic_node<T,T_vector>
  {
  private:
  public:
    logistic_node()
      :stochastic_node<T,T_vector>(2,0)
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      T mu=this->parent(0);
      T sigma=this->parent(1);
      return -(x-mu)/sigma-sigma-2*std::log(1+std::exp(-(x-mu)/sigma));
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      T mu=this->parent(0);
      T sigma=std::abs(this->parent(1));

      return std::make_pair(mu-10*sigma,mu+10*sigma);
      //return make_pair(-10,10);
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,this->parent(0));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new logistic_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };
  
  template <typename T,template <typename TE> class T_vector>
  class logistic_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    logistic_node_factory()
      :abstract_node_factory<T,T_vector>({"mu","sigma"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new logistic_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
}

#endif
