#ifndef LOG_NODE_HPP
#define LOG_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class log_node
    :public cached_dtm_node<T>
  {
  public:
    log_node()
      :cached_dtm_node<T>(1,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return std::log(parent[0]);
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      return std::shared_ptr<node<T> >(new log_node);
    }

  };


  template <typename T>
  class log_vnode
    :public vnode<T>
  {
  public:
    log_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("log",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new log_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new log_vnode<T>(*this));
    }
  };

  template <typename T>
  class log_node_factory
    :public abstract_node_factory<T>
  {
  public:
    log_node_factory()
      :abstract_node_factory<T>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new log_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T>
  log_vnode<T> vlog(const vnode<T>& n1)
  {
    auto result= log_vnode<T>(std::string("log")+node_count<log_vnode<T> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
