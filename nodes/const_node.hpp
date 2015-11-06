#ifndef CONST_NODE_HPP
#define CONST_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class const_node
    :public deterministic_node<T_p,T_var1>
  {
  private:
    T_var1 v;
  public:
    const_node(T_var1 v1)
      :deterministic_node<T_p,T_var1>(0,1),v(v1)
    {
    }
    
    T_var1 do_value(size_t idx)const override
    {
      return v;
    }
  };

  template <typename T_p,typename T_var1>
  class const_vnode
    :public vnode<T_p,T_var1>
  {
  public:
    T_var1 value;
  public:
    const_vnode(std::string n,T_var1 v)
      :vnode<T_p,T_var1>("const",n),value(v)
    {
      this->binded=true;
    }
    
  public:
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new const_node<T_p,T_var1>(value));
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new const_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  class const_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    const_node_factory()
      :abstract_node_factory<T_p,T_var1>({},{"value"},{"value"},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(
	     const std::vector<T_var1>& scalar_param,
	     const std::vector<std::vector<T_var1> >& vector_param)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new const_node<T_p,T_var1>(scalar_param[0]));
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };
  
  using vconst=const_vnode<double,double>;
};

#endif
