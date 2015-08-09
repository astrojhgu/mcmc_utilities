#ifndef PHI_NODE_HPP
#define PHI_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class phi_node
    :public deterministic_node<T_p,T_var1>
  {
  public:
    phi_node()
      :deterministic_node<T_p,T_var1>(1,1)
    {}

    T_var1 do_value(size_t idx,size_t obsid)const override
    {
      return phi(this->parent(0,obsid));
    }
  };


  template <typename T_p,typename T_var1>
  class phi_vnode
    :public vnode<T_p,T_var1>
  {
  public:
    phi_vnode(std::string n,
	      const std::pair<const vnode<T_p,T_var1>&,size_t>& p)
      :vnode<T_p,T_var1>("phi",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new phi_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new phi_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  auto phi(const vnode<T_p,T_var1>& n1)
  {
    return phi_vnode<T_p,T_var1>(std::string("phi")+node_count<phi_vnode<T_p,T_var1> >(),{n1,(size_t)0});
  }
}


#endif
