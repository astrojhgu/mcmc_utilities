#ifndef TAN_NODE_HPP
#define TAN_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class tan_node
    :public deterministic_node<T>
  {
  public:
    tan_node()
      :deterministic_node<T>(1,1)
    {}

    T do_value(size_t idx)const override
    {
      return std::tan(this->parent(0));
    }
  };


  template <typename T>
  class tan_vnode
    :public vnode<T>
  {
  public:
    tan_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("tan",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new tan_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new tan_vnode<T>(*this));
    }
  };

  template <typename T>
  class tan_node_factory
    :public abstract_node_factory<T>
  {
  public:
    tan_node_factory()
      :abstract_node_factory<T>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new tan_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T>
  tan_vnode<T> vtan(const vnode<T>& n1)
  {
    auto result= tan_vnode<T>(std::string("tan")+node_count<tan_vnode<T> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
