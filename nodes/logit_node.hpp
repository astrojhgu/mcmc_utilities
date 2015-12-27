#ifndef LOGIT_NODE_HPP
#define LOGIT_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class logit_node
    :public deterministic_node<T>
  {
  public:
    logit_node()
      :deterministic_node<T>(1,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return logit(parent[0]);
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      return std::shared_ptr<node<T> >(new logit_node);
    }

  };


  template <typename T>
  class logit_vnode
    :public vnode<T>
  {
  public:
    logit_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("logit",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new logit_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new logit_vnode<T>(*this));
    }
  };

  template <typename T>
  class logit_node_factory
    :public abstract_node_factory<T>
  {
  public:
    logit_node_factory()
      :abstract_node_factory<T>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new logit_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T>
  logit_vnode<T> vlogit(const vnode<T>& n1)
  {
    auto result= logit_vnode<T>(std::string("logit")+node_count<logit_vnode<T> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
