#ifndef SQRT_NODE_HPP
#define SQRT_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class sqrt_node
    :public deterministic_node<T_p,T_var1>
  {
  public:
    sqrt_node()
      :deterministic_node<T_p,T_var1>(1,1)
    {}

    T_var1 do_value(size_t idx)const override
    {
      return std::sqrt(this->parent(0));
    }
  };


  template <typename T_p,typename T_var1>
  class sqrt_vnode
    :public vnode<T_p,T_var1>
  {
  public:
    sqrt_vnode(std::string n,
	      const std::pair<const vnode<T_p,T_var1>&,size_t>& p)
      :vnode<T_p,T_var1>("sqrt",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new sqrt_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new sqrt_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  class sqrt_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    sqrt_node_factory()
      :abstract_node_factory<T_p,T_var1>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(const std::vector<T_var1>& hparam)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new sqrt_node<T_p,T_var1>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T_p,typename T_var1>
  sqrt_vnode<T_p,T_var1> vsqrt(const vnode<T_p,T_var1>& n1)
  {
    auto result= sqrt_vnode<T_p,T_var1>(std::string("sqrt")+node_count<sqrt_vnode<T_p,T_var1> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
