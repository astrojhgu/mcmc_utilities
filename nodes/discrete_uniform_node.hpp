#ifndef DISCRETE_UNIFORM_NODE_HPP
#define DISCRETE_UNIFORM_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <helper/abstract_node_factory.hpp>
#include <string>
#include <cmath>

namespace mcmc_utilities
{
  template <typename T>
  class discrete_uniform_node
    :public stochastic_node<T>
  {
  private:
    int a;
    int b;
    std::vector<T> candidates;
  public:
    discrete_uniform_node(int _a,int _b)
      :stochastic_node<T>(0,round((_a+_b)/2)),a(_a),b(_b),candidates(_b-_a)
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

    void do_initialize(size_t n)override
    {
      this->set_value(0,std::round((a+b)/2));
    }
#if 0
    std::vector<T> do_candidate_points()const
    {
      return candidates;
    }
#endif

    T do_regulate(size_t idx,const T& x)const override
    {
      return std::round(x);
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new discrete_uniform_node(a,b);
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T> >(p);
    }

  };

  template <typename T>
  class discrete_uniform_node_factory
    :public abstract_node_factory<T>
  {
  public:
    discrete_uniform_node_factory()
      :abstract_node_factory<T>({},{"x"},{"a","b"})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new discrete_uniform_node<T>(hparam.at(0),hparam.at(1)));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };
  
  template <typename T>
  class discrete_uniform_vnode
    :public vnode<T>
  {
    int a;
    int b;
  public:
    discrete_uniform_vnode(std::string n,int _a,int _b)
      :vnode<T>("discrete_uniform",n),a(_a),b(_b)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new discrete_uniform_node<T>(a,b));
    }
    
    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new discrete_uniform_vnode<T>(*this));
    }
  };

  using vdiscrete_uniform=discrete_uniform_vnode<double>;
};

#endif
