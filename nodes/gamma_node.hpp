#ifndef GAMMA_NODE_HPP
#define GAMMA_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <cmath>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T>
  class gamma_node
    :public stochastic_node<T>
  {
  private:
  public:
    gamma_node()
      :stochastic_node<T>(2,0)
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

    void do_initialize(size_t n) override
    {
      this->set_value(0,this->parent(0)/this->parent(1));
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new gamma_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T> >(p);
    }

  };
  
  
  template <typename T>
  class gamma_vnode
    :public vnode<T>
  {
  public:
    gamma_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("gamma",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new gamma_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new gamma_vnode<T>(*this));
    }
  };

  template <typename T>
  class gamma_node_factory
    :public abstract_node_factory<T>
  {
  public:
    gamma_node_factory()
      :abstract_node_factory<T>({"r","lambda"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new gamma_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vgamma=gamma_vnode<double>;
};

#endif
