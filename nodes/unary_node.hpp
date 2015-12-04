#ifndef UNARY_NODE_HPP
#define UNARY_NODE_HPP

#include <core/deterministic_node.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <functional>


namespace mcmc_utilities
{
  template <typename T>
  class unary_node
    :public deterministic_node<T>
  {
  private:
    std::function<T (const T&)> func;
  public:
    unary_node(const std::function<T (const T&)>& f)
      :deterministic_node<T>(1,1),func(f)
    {}

    T do_value(size_t idx)const override
    {
      return func(this->parent(0));
    }
  };


  template <typename T>
  class unary_vnode
    :public vnode<T>
  {
  private:
    std::function<T (const T&)> func;
  public:
    unary_vnode(std::string n,
		const std::pair<const vnode<T>&,size_t>& p,
		const std::function<T (const T&)>& f)
      :vnode<T>("unary",n,{p}),func(f)
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new unary_node<T>(func));
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new unary_vnode<T>(*this));
    }
  };

  template <typename T>
  unary_vnode<T> vunary(const vnode<T>& n1,const std::function<T (const T&)>& func)
  {
    auto result= unary_vnode<T>(std::string("unary")+node_count<unary_vnode<T> >(),{n1,(size_t)0},func);
    result.named=false;
    return result;
  }

}


#endif
