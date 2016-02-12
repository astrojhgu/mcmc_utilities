#ifndef STEP_NODE_HPP
#define STEP_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class step_node
    :public cached_dtm_node<T>
  {
  public:
    step_node()
      :cached_dtm_node<T>(1,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      //return std::step(parent[0]);
      return static_cast<T>(parent[0]>=0);
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      return std::shared_ptr<node<T> >(new step_node);
    }

  };


  template <typename T>
  class step_vnode
    :public vnode<T>
  {
  public:
    step_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("step",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new step_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new step_vnode<T>(*this));
    }
  };

  template <typename T>
  class step_node_factory
    :public abstract_node_factory<T>
  {
  public:
    step_node_factory()
      :abstract_node_factory<T>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new step_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T>
  step_vnode<T> vstep(const vnode<T>& n1)
  {
    auto result= step_vnode<T>(std::string("step")+node_count<step_vnode<T> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
