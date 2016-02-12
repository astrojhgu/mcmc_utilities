#ifndef EZ_NODE_HPP
#define EZ_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class ez_node
    :public cached_dtm_node<T>
  {
  public:
    ez_node()
      :cached_dtm_node<T>(6,1)
    {}

    T E(T z,T Omega_m,T Omega_l,T Omega_k,T Omega_rad,T w)const
    {
      T zp1=1+z;
      return std::sqrt(Omega_rad*zp1*zp1*zp1*zp1+Omega_m*zp1*zp1*zp1+Omega_k*zp1*zp1+Omega_l*pow(zp1,3*(1+w)));
    }

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      T z=parent[0];
      T Omega_m=parent[1];
      T Omega_l=parent[2];
      T Omega_k=parent[3];
      T Omega_rad=parent[4];
      T w=parent[5];
      return E(z,Omega_m,Omega_l,Omega_k,Omega_rad,w);
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      return std::shared_ptr<node<T> >(new ez_node);
    }

  };


  template <typename T>
  class ez_vnode
    :public vnode<T>
  {
  public:
    ez_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("ez",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new ez_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new ez_vnode<T>(*this));
    }
  };

  template <typename T>
  class ez_node_factory
    :public abstract_node_factory<T>
  {
  public:
    ez_node_factory()
      :abstract_node_factory<T>({"z","Omega_m","Omega_l","Omega_k","Omega_rad","w"},{"E"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new ez_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T>
  ez_vnode<T> vez(const vnode<T>& n1)
  {
    auto result= ez_vnode<T>(std::string("ez")+node_count<ez_vnode<T> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
