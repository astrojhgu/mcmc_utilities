#ifndef ILOGIT_NODE_HPP
#define ILOGIT_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class ilogit_node
    :public deterministic_node<T>
  {
  public:
    ilogit_node()
      :deterministic_node<T>(1,1)
    {}

    T do_value(size_t idx)const override
    {
      return ilogit(this->parent(0));
    }
  };


  template <typename T>
  class ilogit_vnode
    :public vnode<T>
  {
  public:
    ilogit_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("ilogit",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new ilogit_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new ilogit_vnode<T>(*this));
    }
  };

  template <typename T>
  class ilogit_node_factory
    :public abstract_node_factory<T>
  {
  public:
    ilogit_node_factory()
      :abstract_node_factory<T>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new ilogit_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T>
  ilogit_vnode<T> vilogit(const vnode<T>& n1)
  {
    auto result= ilogit_vnode<T>(std::string("ilogit")+node_count<ilogit_vnode<T> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
