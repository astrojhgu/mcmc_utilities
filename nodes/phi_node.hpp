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
  template <typename T>
  class phi_node
    :public deterministic_node<T>
  {
  public:
    phi_node()
      :deterministic_node<T>(1,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return phi(parent[0]);
    }
  };


  template <typename T>
  class phi_vnode
    :public vnode<T>
  {
  public:
    phi_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("phi",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new phi_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new phi_vnode<T>(*this));
    }
  };

  template <typename T>
  class phi_node_factory
    :public abstract_node_factory<T>
  {
  public:
    phi_node_factory()
      :abstract_node_factory<T>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new phi_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T>
  phi_vnode<T> vphi(const vnode<T>& n1)
  {
    auto result= phi_vnode<T>(std::string("phi")+node_count<phi_vnode<T> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
