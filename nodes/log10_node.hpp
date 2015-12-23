#ifndef LOG10_NODE_HPP
#define LOG10_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class log10_node
    :public deterministic_node<T>
  {
  public:
    log10_node()
      :deterministic_node<T>(1,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return std::log10(parent[0]);
    }
  };


  template <typename T>
  class log10_vnode
    :public vnode<T>
  {
  public:
    log10_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("log10",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new log10_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new log10_vnode<T>(*this));
    }
  };

  template <typename T>
  class log10_node_factory
    :public abstract_node_factory<T>
  {
  public:
    log10_node_factory()
      :abstract_node_factory<T>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new log10_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T>
  log10_vnode<T> vlog10(const vnode<T>& n1)
  {
    auto result= log10_vnode<T>(std::string("log10")+node_count<log10_vnode<T> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
