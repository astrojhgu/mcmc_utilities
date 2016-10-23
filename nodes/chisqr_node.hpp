#ifndef CHISQR_NODE_HPP
#define CHISQR_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/node_counter.hpp>
#include <limits>
#include <string>
#include <cmath>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class chisqr_node
    :public stochastic_node<T,T_vector>
  {
  private:
    int k;
  public:
    chisqr_node(int k1)
      :stochastic_node<T,T_vector>(0,static_cast<T>(k)),
      k(k1)
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      return (static_cast<T>(k)/2-1)*std::log(x)-x/2-(static_cast<T>(k)/2)*std::log(2)-lgamma(static_cast<T>(k)/2);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      return std::make_pair(std::numeric_limits<T>::epsilon(),T(10.0*k));
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,static_cast<T>(k));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new chisqr_node(k);
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class chisqr_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    chisqr_node_factory()
      :abstract_node_factory<T,T_vector>({},{"x"},{"k"})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new chisqr_node<T,T_vector>(hparam[0]));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };
}

#endif
