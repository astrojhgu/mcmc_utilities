#ifndef NORMAL_NODE_HPP
#define NORMAL_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class normal_node
    :public stochastic_node<T_p,T_var1>
  {
  private:
  public:
    normal_node()
      :stochastic_node<T_p,T_var1>(2,0)
    {}
    
  private:
    T_p do_log_prob()const override
    {
      const static T_var1 PI=std::atan(1.0)*4;
      T_var1 x=this->value(0);
      T_var1 mu=this->parent(0);
      T_var1 sigma=this->parent(1);
      return -(x-mu)*(x-mu)/(2*sigma*sigma)-std::log(sigma*std::sqrt(2*PI));
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T_var1,T_var1> do_var_range()const override
    {
      T_var1 mu=this->parent(0);
      T_var1 sigma=std::abs(this->parent(1));

      return std::make_pair(mu-5*sigma,mu+5*sigma);
      //return make_pair(-10,10);
    }
  };
  
  
  template <typename T_p,typename T_var1>
  class normal_vnode
    :public vnode<T_p,T_var1>
  {
  public:
    normal_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T_p,T_var1>&,size_t> >& p)
      :vnode<T_p,T_var1>("normal",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new normal_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new normal_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  class normal_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    normal_node_factory()
      :abstract_node_factory<T_p,T_var1>({"mu","sigma"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(const std::vector<T_var1>& hparam)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new normal_node<T_p,T_var1>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vnormal=normal_vnode<double,double>;
};

#endif
