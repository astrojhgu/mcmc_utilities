#ifndef EZ_NODE_HPP
#define EZ_NODE_HPP

#include <core/differentiable_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class ez_node
    :public differentiable_dtm_node<T,T_vector>
  {
  public:
    ez_node()
      :differentiable_dtm_node<T,T_vector>(6,1)
    {}

    T E(T z,T Omega_m,T Omega_l,T Omega_k,T Omega_rad,T w)const
    {
      T zp1=1+z;
      return std::sqrt(Omega_rad*zp1*zp1*zp1*zp1+Omega_m*zp1*zp1*zp1+Omega_k*zp1*zp1+Omega_l*pow(zp1,3*(1+w)));
    }

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      T z=parent[0];
      T Omega_m=parent[1];
      T Omega_l=parent[2];
      T Omega_k=parent[3];
      T Omega_rad=parent[4];
      T w=parent[5];
      return E(z,Omega_m,Omega_l,Omega_k,Omega_rad,w);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new ez_node);
    }
  };


  template <typename T,template <typename TE> class T_vector>
  class ez_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    ez_node_factory()
      :abstract_node_factory<T,T_vector>({"z","Omega_m","Omega_l","Omega_k","Omega_rad","w"},{"E"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new ez_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };
}


#endif
