#ifndef NORMAL_NODE_HPP
#define NORMAL_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T>
  class normal_node
    :public stochastic_node<T>
  {
  private:
  public:
    normal_node()
      :stochastic_node<T>(2,0)
    {}
    
  private:
    T do_log_prob()const override
    {
      const static T PI=std::atan(1.0)*4;
      T x=this->value(0);
      T mu=this->parent(0);
      T sigma=this->parent(1);
      return -(x-mu)*(x-mu)/(2*sigma*sigma)-std::log(sigma*std::sqrt(2*PI));
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      T mu=this->parent(0);
      T sigma=std::abs(this->parent(1));

      return std::make_pair(mu-5*sigma,mu+5*sigma);
      //return make_pair(-10,10);
    }
  };
  
  
  template <typename T>
  class normal_vnode
    :public vnode<T>
  {
  public:
    normal_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("normal",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new normal_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new normal_vnode<T>(*this));
    }
  };

  template <typename T>
  class normal_node_factory
    :public abstract_node_factory<T>
  {
  public:
    normal_node_factory()
      :abstract_node_factory<T>({"mu","sigma"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new normal_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vnormal=normal_vnode<double>;
};

#endif
