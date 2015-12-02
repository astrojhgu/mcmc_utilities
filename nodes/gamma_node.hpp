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
  template <typename T_p,typename T_var1>
  class gamma_node
    :public stochastic_node<T_p,T_var1>
  {
  private:
  public:
    gamma_node()
      :stochastic_node<T_p,T_var1>(2,0)
    {}
    
  private:
    T_p do_log_prob()const override
    {
      T_var1 x=this->value(0);
      T_var1 r=this->parent(0);
      T_var1 lambda=this->parent(1);
      return r*std::log(lambda)+(r-1)*std::log(x)-lambda*x-std::lgamma(r);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T_var1,T_var1> do_var_range()const override
    {
      return std::make_pair((T_var1)0,(T_var1)(20.)/this->parent(1));
    }
  };
  
  
  template <typename T_p,typename T_var1>
  class gamma_vnode
    :public vnode<T_p,T_var1>
  {
  public:
    gamma_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T_p,T_var1>&,size_t> >& p)
      :vnode<T_p,T_var1>("gamma",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new gamma_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new gamma_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  class gamma_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    gamma_node_factory()
      :abstract_node_factory<T_p,T_var1>({"r","lambda"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(const std::vector<T_var1>& hparam)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new gamma_node<T_p,T_var1>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vgamma=gamma_vnode<double,double>;
};

#endif
