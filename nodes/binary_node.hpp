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
  template <typename T_p,typename T_var1>
  class binary_node
    :public deterministic_node<T_p,T_var1>
  {
  private:
    std::function<T_var1 (const T_var1&)> func;
  public:
    binary_node(const std::function<T_var1 (const T_var1&,const T_var1&)>& f)
      :deterministic_node<T_p,T_var1>(1,1),func(f)
    {}

    T_var1 do_value(size_t idx)const override
    {
      return func(this->parent(0));
    }
  };


  template <typename T_p,typename T_var1>
  class binary_vnode
    :public vnode<T_p,T_var1>
  {
  private:
    std::function<T_var1 (const T_var1&,const T_var1&)> func;
  public:
    binary_vnode(std::string n,
		 const std::pair<const vnode<T_p,T_var1>&,size_t>& p,
		 const std::function<T_var1 (const T_var1&,const T_var1&)>& f)
      :vnode<T_p,T_var1>("binary",n,{p}),func(f)
    {
      this->binded=true;
    }

    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new binary_node<T_p,T_var1>(func));
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new binary_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  binary_vnode<T_p,T_var1> vbinary(const vnode<T_p,T_var1>& n1,const std::function<T_var1 (const T_var1&,const T_var1&)>& func)
  {
    auto result= binary_vnode<T_p,T_var1>(std::string("binary")+node_count<binary_vnode<T_p,T_var1> >(),{n1,(size_t)0},func);
    result.named=false;
    return result;
  }

}


#endif
