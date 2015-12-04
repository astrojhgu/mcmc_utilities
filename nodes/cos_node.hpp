#ifndef COS_NODE_HPP
#define COS_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class cos_node
    :public deterministic_node<T_p,T_var1>
  {
  public:
    cos_node()
      :deterministic_node<T_p,T_var1>(1,1)
    {}

    T_var1 do_value(size_t idx)const override
    {
      return std::cos(this->parent(0));
    }
  };


  template <typename T_p,typename T_var1>
  class cos_vnode
    :public vnode<T_p,T_var1>
  {
  public:
    cos_vnode(std::string n,
	      const std::pair<const vnode<T_p,T_var1>&,size_t>& p)
      :vnode<T_p,T_var1>("cos",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new cos_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new cos_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  class cos_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    cos_node_factory()
      :abstract_node_factory<T_p,T_var1>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(const std::vector<T_var1>& hparam)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new cos_node<T_p,T_var1>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T_p,typename T_var1>
  cos_vnode<T_p,T_var1> vcos(const vnode<T_p,T_var1>& n1)
  {
    auto result= cos_vnode<T_p,T_var1>(std::string("cos")+node_count<cos_vnode<T_p,T_var1> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
