#ifndef BETA_NODE_HPP
#define BETA_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <limits>
#include <string>
#include <cmath>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T>
  class beta_node
    :public stochastic_node<T>
  {
  private:
    T alpha,beta;
  public:
    beta_node(T a,T b)
      :stochastic_node<T>(0,.5),alpha(a),beta(b)
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      return (alpha-1)*std::log(x)+(beta-1)*std::log(1-x)-lbeta(alpha,beta);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      return std::make_pair(std::numeric_limits<T>::epsilon(),T(1-std::numeric_limits<T>::epsilon()));
    }
  };
  
  
  template <typename T>
  class beta_vnode
    :public vnode<T>
  {
  private:
    T a,b;
  public:
    beta_vnode(std::string n,const T& v1,const T& v2)
      :vnode<T>("beta",n),a(v1),b(v2)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new beta_node<T>(a,b));
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new beta_vnode<T>(*this));
    }
  };

  template <typename T>
  class beta_node_factory
    :public abstract_node_factory<T>
  {
  public:
    beta_node_factory()
      :abstract_node_factory<T>({},{"x"},{"a","b"})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new beta_node<T>(hparam[0],hparam[1]));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vbeta=beta_vnode<double>;
};

#endif
