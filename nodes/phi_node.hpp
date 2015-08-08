#ifndef PHI_NODE_HPP
#define PHI_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <initializer_list>


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
  class _phi_vnode
    :public _vnode<T_p,T_var1>
  {
  public:
    _phi_vnode(std::string n,const std::initializer_list<std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> >& p)
      :_vnode<T_p,T_var1>("phi",n,p)
    {
      this->binded=true;
    }

    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new phi_node<T_p,T_var1>);
    }
  };

  template <typename T_p,typename T_var1>
  auto phi(const cpnt<T_p,T_var1>& n1)
  {
    return cpnt<T_p,T_var1>(std::shared_ptr<_vnode<T_p,T_var1> >(
								 new _phi_vnode<T_p,T_var1>(std::string("phi")+node_count<_add_vnode<T_p,T_var1> >(),{{n1.pn,n1.n}})));
  }
}


#endif
