#ifndef PARETO_NODE_HPP
#define PARETO_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>
#include <math/distributions.hpp>

namespace mcmc_utilities
{
  template <typename T>
  class pareto_node
    :public stochastic_node<T>
  {
  private:
  public:
    pareto_node()
      :stochastic_node<T>(2,0)
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
      
      return std::make_pair(c,std::exp(6/(alpha+1));
    }

    void do_initialize(size_t n)override
    {
      this->set_value(0,this->parent(0));
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new pareto_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T> >(p);
    }

  };
  
  
  template <typename T>
  class pareto_vnode
    :public vnode<T>
  {
  public:
    pareto_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("pareto",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new pareto_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new pareto_vnode<T>(*this));
    }
  };

  template <typename T>
  class pareto_node_factory
    :public abstract_node_factory<T>
  {
  public:
    pareto_node_factory()
      :abstract_node_factory<T>({"mu","sigma"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new pareto_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vpareto=pareto_vnode<double>;
};

#endif
