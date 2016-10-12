#ifndef CONST_NODE_HPP
#define CONST_NODE_HPP
#include <core/deterministic_node.hpp>
#include <core/tp_aware_dtm_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <helper/abstract_node_factory.hpp>
#include <string>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class const_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  private:
    T v;
  public:
    const_node(T v1)
      :tp_aware_dtm_node<T,T_vector>(0,1),v(v1)
    {
    }
    
    T do_calc(size_t idx,const T_vector<T>&)const override
    {
      return v;
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const
    {
      return order{0,true,true};
    }
    
    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new const_node(v));
    }

    void set_value(T v1)
    {
      v=v1;
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class const_vnode
    :public vnode<T,T_vector>
  {
  public:
    T value;
  public:
    const_vnode(std::string n,T v)
      :vnode<T,T_vector>("const",n),value(v)
    {
      this->binded=true;
    }
    
  public:
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new const_node<T,T_vector>(value));
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new const_vnode<T,T_vector>(*this));
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class const_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    const_node_factory()
      :abstract_node_factory<T,T_vector>({},{"v"},{"value"})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new const_node<T,T_vector>(hparam[0]));
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };
  
  //using vconst=const_vnode<double>;
};

#endif
