#ifndef SIN_NODE_HPP
#define SIN_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class sin_node
    :public deterministic_node<T>
  {
  public:
    sin_node()
      :deterministic_node<T>(1,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return std::sin(parent[0]);
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new sin_node;
      return std::shared_ptr<node<T> >(p);
    }

  };


  template <typename T>
  class sin_vnode
    :public vnode<T>
  {
  public:
    sin_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("sin",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new sin_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new sin_vnode<T>(*this));
    }
  };

  template <typename T>
  class sin_node_factory
    :public abstract_node_factory<T>
  {
  public:
    sin_node_factory()
      :abstract_node_factory<T>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new sin_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T>
  sin_vnode<T> vsin(const vnode<T>& n1)
  {
    auto result= sin_vnode<T>(std::string("sin")+node_count<sin_vnode<T> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
