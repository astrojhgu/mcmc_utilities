#ifndef COSMOLOGY_NODE_HPP
#define COSMOLOGY_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <coscalc.hpp>


namespace mcmc_utilities
{

  ///////////luminosity distance node//////////////
  template <typename T,template <typename TE> class T_vector>
  class luminosity_distance_node
    :public cached_dtm_node<T,T_vector>
  {
  private:
    struct parent_value_cache_t
    {
    public:
      T values[7];
      T& z;
      T& H0;
      T& Omega_m;
      T& Omega_l;
      T& Omega_k;
      T& Omega_rad;
      T& w;
      
      parent_value_cache_t()
	:z(values[0]),H0(values[1]),Omega_m(values[2]),
	 Omega_l(values[3]),Omega_k(values[4]),Omega_rad(values[5]),w(values[6])
      {}
    }parent_value_cache;
    T value_cache;
  public:
    luminosity_distance_node()
      :cached_dtm_node<T,T_vector>(7,1),parent_value_cache{},value_cache{}
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
#if 0
      bool changed=false;
      for(int i=0;i<6;++i)
	{
	  T x=parent[i];
	  if(parent_value_cache.values[i]!=x)
	    {
	      const_cast<parent_value_cache_t&>(parent_value_cache).values[i]=x;
	      changed=true;
	    }
	}

      if(!changed)
	{
	  return value_cache;
	}
      
      
      T z=parent_value_cache.z;
      T H0=parent_value_cache.H0;
      T Omega_m=parent_value_cache.Omega_m;
      T Omega_l=parent_value_cache.Omega_l;
      T Omega_k=parent_value_cache.Omega_k;
      T Omega_rad=parent_value_cache.Omega_rad;
      T w=parent_value_cache.w;
#else
      T z=parent[0];
      T H0=parent[1];
      T Omega_m=parent[2];
      T Omega_l=parent[3];
      T Omega_k=parent[4];
      T Omega_rad=parent[5];
      T w=parent[6];
#endif
      
      
      using u=coscalc::unit<coscalc::unit_system::SI>;
      coscalc::cosmic_param cp(2.99792458E8*u::m/u::s,H0*u::km/u::s/u::Mpc,
			       Omega_m,Omega_l,Omega_k,Omega_rad,w);
      coscalc::cosmic_calculator cc(cp);
      //const_cast<T&>(value_cache)=cc.calc_dl(z);
      //return value_cache;
      return cc.calc_dl(z);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new luminosity_distance_node);
    }

  };


  template <typename T,template <typename TE> class T_vector>
  class luminosity_distance_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    luminosity_distance_node_factory()
      :abstract_node_factory<T,T_vector>({"z","H0","Omega_m","Omega_l","Omega_k","Omega_rad","w"},{"D_L"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new luminosity_distance_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  ///////////asize distance node//////////////
  template <typename T,template <typename TE> class T_vector>
  class asize_distance_node
    :public cached_dtm_node<T,T_vector>
  {
  private:
    struct parent_value_cache_t
    {
    public:
      T values[7];
      T& z;
      T& H0;
      T& Omega_m;
      T& Omega_l;
      T& Omega_k;
      T& Omega_rad;
      T& w;
      
      parent_value_cache_t()
	:z(values[0]),H0(values[1]),Omega_m(values[2]),
	 Omega_l(values[3]),Omega_k(values[4]),Omega_rad(values[5]),w(values[6])
      {}
    }parent_value_cache;
    T value_cache;

  public:
    asize_distance_node()
      :cached_dtm_node<T,T_vector>(7,1),parent_value_cache{},value_cache{}
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      bool changed=false;
      for(int i=0;i<6;++i)
	{
	  T x=parent[i];
	  if(parent_value_cache.values[i]!=x)
	    {
	      const_cast<parent_value_cache_t&>(parent_value_cache).values[i]=x;
	      changed=true;
	    }
	}

      if(!changed)
	{
	  return value_cache;
	}
      
      
      T z=parent_value_cache.z;
      T H0=parent_value_cache.H0;
      T Omega_m=parent_value_cache.Omega_m;
      T Omega_l=parent_value_cache.Omega_l;
      T Omega_k=parent_value_cache.Omega_k;
      T Omega_rad=parent_value_cache.Omega_rad;
      T w=parent_value_cache.w;

      
      
      using u=coscalc::unit<coscalc::unit_system::SI>;
      coscalc::cosmic_param cp(2.99792458E8*u::m/u::s,H0*u::km/u::s/u::Mpc,
			       Omega_m,Omega_l,Omega_k,Omega_rad,w);
      coscalc::cosmic_calculator cc(cp);
      const_cast<T&>(value_cache)=cc.calc_da(z);
      return value_cache;
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new asize_distance_node);
    }

  };

  template <typename T,template <typename TE> class T_vector>
  class asize_distance_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    asize_distance_node_factory()
      :abstract_node_factory<T,T_vector>({"z","H0","Omega_m","Omega_l","Omega_k","Omega_rad","w"},{"D_L"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new asize_distance_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

}


#endif
