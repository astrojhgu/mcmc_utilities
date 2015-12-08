#ifndef COSMOLOGY_NODE_HPP
#define COSMOLOGY_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <coscalc.hpp>


namespace mcmc_utilities
{

  ///////////luminosity distance node//////////////
  template <typename T>
  class luminosity_distance_node
    :public deterministic_node<T>
  {
  public:
    luminosity_distance_node()
      :deterministic_node<T>(5,1)
    {}

    T do_value(size_t idx)const override
    {
      //return std::log(this->parent(0));
      //cosmic_param cp(2.
      T z=this->parent(0);
      T H0=this->parent(1);
      T Omega_m=this->parent(2);
      T Omega_l=this->parent(3);
      T Omega_k=this->parent(4);
      using u=coscalc::unit<coscalc::unit_system::SI>;
      coscalc::cosmic_param cp(2.99792458E8*u::m/u::s,H0*u::km/u::s/u::Mpc,
		      Omega_m,Omega_l,Omega_k);
      coscalc::cosmic_calculator cc(cp);
      return cc.calc_dl(z);
    }
  };


  template <typename T>
  class luminosity_distance_vnode
    :public vnode<T>
  {
  public:
    luminosity_distance_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("luminosity_distance",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new luminosity_distance_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new luminosity_distance_vnode<T>(*this));
    }
  };

  template <typename T>
  class luminosity_distance_node_factory
    :public abstract_node_factory<T>
  {
  public:
    luminosity_distance_node_factory()
      :abstract_node_factory<T>({"z","H0","Omega_m","Omega_l","Omega_k"},{"D_L"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new luminosity_distance_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  ///////////asize distance node//////////////
  template <typename T>
  class asize_distance_node
    :public deterministic_node<T>
  {
  public:
    asize_distance_node()
      :deterministic_node<T>(5,1)
    {}

    T do_value(size_t idx)const override
    {
      //return std::log(this->parent(0));
      //cosmic_param cp(2.
      T z=this->parent(0);
      T H0=this->parent(1);
      T Omega_m=this->parent(2);
      T Omega_l=this->parent(3);
      T Omega_k=this->parent(4);
      using u=coscalc::unit<coscalc::unit_system::SI>;
      coscalc::cosmic_param cp(2.99792458E8*u::m/u::s,H0*u::km/u::s/u::Mpc,
		      Omega_m,Omega_l,Omega_k);
      coscalc::cosmic_calculator cc(cp);
      return cc.calc_da(z);
    }
  };


  template <typename T>
  class asize_distance_vnode
    :public vnode<T>
  {
  public:
    asize_distance_vnode(std::string n,
	      const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("asize_distance",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new asize_distance_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new asize_distance_vnode<T>(*this));
    }
  };

  template <typename T>
  class asize_distance_node_factory
    :public abstract_node_factory<T>
  {
  public:
    asize_distance_node_factory()
      :abstract_node_factory<T>({"z","H0","Omega_m","Omega_l","Omega_k"},{"D_L"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new asize_distance_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

}


#endif
