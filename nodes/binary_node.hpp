#ifndef BINARY_NODE_HPP
#define BINARY_NODE_HPP

#include <core/deterministic_node.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <functional>


namespace mcmc_utilities
{
  template <typename T>
  class binary_node
    :public deterministic_node<T>
  {
  private:
    std::function<T (const T&)> func;
  public:
    binary_node(const std::function<T (const T&,const T&)>& f)
      :deterministic_node<T>(1,1),func(f)
    {}

    T do_value(size_t idx)const override
    {
      return func(this->parent(0));
    }
  };


  template <typename T>
  class binary_vnode
    :public vnode<T>
  {
  private:
    std::function<T (const T&,const T&)> func;
  public:
    binary_vnode(std::string n,
		 const std::pair<const vnode<T>&,size_t>& p,
		 const std::function<T (const T&,const T&)>& f)
      :vnode<T>("binary",n,{p}),func(f)
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new binary_node<T>(func));
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new binary_vnode<T>(*this));
    }
  };

  template <typename T>
  binary_vnode<T> vbinary(const vnode<T>& n1,const std::function<T (const T&,const T&)>& func)
  {
    auto result= binary_vnode<T>(std::string("binary")+node_count<binary_vnode<T> >(),{n1,(size_t)0},func);
    result.named=false;
    return result;
  }

}


#endif
