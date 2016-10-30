#ifndef EXP_DIST_NODE_HPP
#define EXP_DIST_NODE_HPP

//double exponential distribution

#include <core/forward_sampleable_node.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class exp_dist_node
    :public forward_sampleable_node<T,T_vector>
  {
  private:
  public:
    exp_dist_node()
      :forward_sampleable_node<T,T_vector>(1,std::numeric_limits<T>::epsilon())
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      T lambda=this->parent(0);
      return std::log(lambda)-lambda*x;
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      T lambda=this->parent(0);
      return std::make_pair(0,static_cast<T>(20)/lambda);
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,static_cast<T>(1)/this->parent(0));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new exp_dist_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  class exp_dist_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    exp_dist_node_factory()
      :abstract_node_factory<T,T_vector>({"lambda"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new exp_dist_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };
}

#endif
