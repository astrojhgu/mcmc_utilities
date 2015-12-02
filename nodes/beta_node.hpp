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
  template <typename T_p,typename T_var1>
  class beta_node
    :public stochastic_node<T_p,T_var1>
  {
  private:
    T_var1 alpha,beta;
  public:
    beta_node(T_var1 a,T_var1 b)
      :stochastic_node<T_p,T_var1>(0,.5),alpha(a),beta(b)
    {}
    
  private:
    T_p do_log_prob()const override
    {
      T_var1 x=this->value(0);
      return (alpha-1)*std::log(x)+(beta-1)*std::log(1-x)-lbeta(alpha,beta);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T_var1,T_var1> do_var_range()const override
    {
      return std::make_pair(std::numeric_limits<T_var1>::epsilon(),T_var1(1-std::numeric_limits<T_var1>::epsilon()));
    }
  };
  
  
  template <typename T_p,typename T_var1>
  class beta_vnode
    :public vnode<T_p,T_var1>
  {
  private:
    T_var1 a,b;
  public:
    beta_vnode(std::string n,const T_var1& v1,const T_var1& v2)
      :vnode<T_p,T_var1>("beta",n),a(v1),b(v2)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new beta_node<T_p,T_var1>(a,b));
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new beta_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  class beta_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    beta_node_factory()
      :abstract_node_factory<T_p,T_var1>({},{"x"},{"a","b"})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(const std::vector<T_var1>& hparam)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new beta_node<T_p,T_var1>(hparam[0],hparam[1]));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vbeta=beta_vnode<double,double>;
};

#endif
