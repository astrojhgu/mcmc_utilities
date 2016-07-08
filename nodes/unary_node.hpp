#ifndef UNARY_NODE_HPP
#define UNARY_NODE_HPP

#include <core/tp_aware_dtm_node.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <functional>


namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class unary_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  private:
    std::function<T (const T&)> func;
  public:
    unary_node(const std::function<T (const T&)>& f)
      :tp_aware_dtm_node<T,T_vector>(1,1),func(f)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return func(parent[0]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new unary_node(func));
    }

  };


  template <typename T,template <typename TE> class T_vector>
  class unary_vnode
    :public vnode<T,T_vector>
  {
  private:
    std::function<T (const T&)> func;
  public:
    unary_vnode(std::string n,
		const std::pair<const vnode<T,T_vector>&,size_t>& p,
		const std::function<T (const T&)>& f)
      :vnode<T,T_vector>("unary",n,{p}),func(f)
    {
      this->binded=true;
    }

    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new unary_node<T,T_vector>(func));
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new unary_vnode<T,T_vector>(*this));
    }
  };

  template <typename T,template <typename TE> class T_vector>
  unary_vnode<T,T_vector> vunary(const vnode<T,T_vector>& n1,const std::function<T (const T&)>& func)
  {
    auto result= unary_vnode<T,T_vector>(std::string("unary")+node_count<unary_vnode<T,T_vector> >(),{n1,(size_t)0},func);
    result.named=false;
    return result;
  }

}


#endif
