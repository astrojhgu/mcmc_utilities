#ifndef CONST_NODE_HPP
#define CONST_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <helper/abstract_node_factory.hpp>
#include <string>

namespace mcmc_utilities
{
  template <typename T>
  class const_node
    :public deterministic_node<T>
  {
  private:
    T v;
  public:
    const_node(T v1)
      :deterministic_node<T>(0,1),v(v1)
    {
    }
    
    T do_calc(size_t idx,const std::vector<T>&)const override
    {
      return v;
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      return std::shared_ptr<node<T> >(new const_node(v));
    }

    void set_value(T v1)
    {
      v=v1;
    }
  };

  template <typename T>
  class const_vnode
    :public vnode<T>
  {
  public:
    T value;
  public:
    const_vnode(std::string n,T v)
      :vnode<T>("const",n),value(v)
    {
      this->binded=true;
    }
    
  public:
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new const_node<T>(value));
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new const_vnode<T>(*this));
    }
  };

  template <typename T>
  class const_node_factory
    :public abstract_node_factory<T>
  {
  public:
    const_node_factory()
      :abstract_node_factory<T>({},{"v"},{"value"})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new const_node<T>(hparam[0]));
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };
  
  using vconst=const_vnode<double>;
};

#endif
