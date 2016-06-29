#ifndef LOGIT_NODE_HPP
#define LOGIT_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class logit_node
    :public cached_dtm_node<T,T_vector>
  {
  public:
    logit_node()
      :cached_dtm_node<T,T_vector>(1,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return logit(parent[0]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new logit_node);
    }

  };


  template <typename T,template <typename TE> class T_vector>
  class logit_vnode
    :public vnode<T,T_vector>
  {
  public:
    logit_vnode(std::string n,
	      const std::pair<const vnode<T,T_vector>&,size_t>& p)
      :vnode<T,T_vector>("logit",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new logit_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new logit_vnode<T,T_vector>(*this));
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class logit_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    logit_node_factory()
      :abstract_node_factory<T,T_vector>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new logit_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T,template <typename TE> class T_vector>
  logit_vnode<T,T_vector> vlogit(const vnode<T,T_vector>& n1)
  {
    auto result= logit_vnode<T,T_vector>(std::string("logit")+node_count<logit_vnode<T,T_vector> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
