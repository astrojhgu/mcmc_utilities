#ifndef DISCRETE_UNIFORM_NODE_HPP
#define DISCRETE_UNIFORM_NODE_HPP
#include <core/forward_sampleable_node.hpp>
#include <helper/node_counter.hpp>
#include <helper/abstract_node_factory.hpp>
#include <string>
#include <cmath>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class discrete_uniform_node
    :public forward_sampleable_node<T,T_vector>
  {
  private:
    int a;
    int b;
    T_vector<T> candidates;
  public:
    discrete_uniform_node(int _a,int _b)
      :forward_sampleable_node<T,T_vector>(0,round((_a+_b)/2)),a(_a),b(_b),candidates(_b-_a)
    {
      for(int i=a;i<b;++i)
	{
	  candidates.at(i-a)=i;
	}
    }
    
  private:
    T do_log_prob()const override
    {
      return 0;
    }
    
    bool is_continuous(size_t)const override
    {
      return false;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      return std::make_pair<T,T>(a,b);
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,std::round((a+b)/2));
    }
#if 0
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
      auto p=new discrete_uniform_node(a,b);
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };

  template <typename T,template <typename TE> class T_vector>
  class discrete_uniform_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    discrete_uniform_node_factory()
      :abstract_node_factory<T,T_vector>({},{"x"},{"a","b"})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new discrete_uniform_node<T,T_vector>(hparam.at(0),hparam.at(1)));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };
};

#endif
