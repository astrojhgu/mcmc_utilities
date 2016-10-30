#ifndef BERN_NODE_HPP
#define BERN_NODE_HPP
#include <core/forward_sampleable_node.hpp>
#include <helper/node_counter.hpp>
#include <helper/abstract_node_factory.hpp>
#include <string>
#include <cmath>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class bern_node
    :public forward_sampleable_node<T,T_vector>
  {
  private:
    T_vector<T> candidates;
  public:
    bern_node()
      :forward_sampleable_node<T,T_vector>(1,0),candidates{0,1}
    {
    }
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      T p=this->parent(0);
      //return std::log(int(x)==0?(1-p):p);
      return std::round(x)==0?std::log(1-p):std::log(p);
    }
    
    bool is_continuous(size_t)const override
    {
      return false;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      return std::make_pair<T,T>(0,2);
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,0);
    }
#if 1
    T_vector<T> do_candidate_points()const
    {
      return candidates;
    }
#endif

    T do_regulate(size_t idx,const T& x)const override
    {
      return std::round(x);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new bern_node();
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };

  template <typename T,template <typename TE> class T_vector>
  class bern_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    bern_node_factory()
      :abstract_node_factory<T,T_vector>({"p"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new bern_node<T,T_vector>());
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };
}

#endif
