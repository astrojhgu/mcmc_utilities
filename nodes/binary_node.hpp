#ifndef BINARY_NODE_HPP
#define BINARY_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <functional>


namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class binary_node
    :public cached_dtm_node<T,T_vector>
  {
  private:
    std::function<T (const T&)> func;
  public:
    binary_node(const std::function<T (const T&,const T&)>& f)
      :cached_dtm_node<T,T_vector>(1,1),func(f)
    {}

    T do_calc(size_t idx,const std::vector<T,T_vector>& parent)const override
    {
      return func(this->parent[0]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new binary_node(alpha,beta);
      return std::shared_ptr<node<T,T_vector> >(p);
    }
    
  };


  template <typename T,template <typename TE> class T_vector>
  class binary_vnode
    :public vnode<T,T_vector>
  {
  private:
    std::function<T (const T&,const T&)> func;
  public:
    binary_vnode(std::string n,
		 const std::pair<const vnode<T,T_vector>&,size_t>& p,
		 const std::function<T (const T&,const T&)>& f)
      :vnode<T,T_vector>("binary",n,{p}),func(f)
    {
      this->binded=true;
    }

    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new binary_node<T,T_vector>(func));
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new binary_vnode<T,T_vector>(*this));
    }
  };

  template <typename T,template <typename TE> class T_vector>
  binary_vnode<T,T_vector> vbinary(const vnode<T,T_vector>& n1,const std::function<T (const T&,const T&)>& func)
  {
    auto result= binary_vnode<T,T_vector>(std::string("binary")+node_count<binary_vnode<T,T_vector> >(),{n1,(size_t)0},func);
    result.named=false;
    return result;
  }

}


#endif
