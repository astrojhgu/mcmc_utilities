#ifndef VNODE_HPP
#define VNODE_HPP
#include "node_counter.hpp"
#include <core/deterministic_node.hpp>
#include <string>
#include <vector>
#include <initializer_list>
#include <map>
#include <memory>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class eq_vnode;

  template <typename T_p,typename T_var1>
  class eq_node;
  
  template <typename T_p,typename T_var1>
  class vnode
  {
  public:
    std::string type;
    std::string name;
    std::vector<std::pair<std::shared_ptr<vnode<T_p,T_var1> >,size_t> > parents;
    bool binded;
    bool added;
    
    vnode(const std::string& n)
      :vnode("unbinded",n,{})
    {}
    
    vnode(const std::string& t,const std::string& n)
      :vnode(t,n,{})
    {}
    
    vnode(const std::string& t,const std::string& n,const std::initializer_list<std::pair<const vnode<T_p,T_var1>&,size_t> >& p)
      :type(t),name(n),binded(false),added(false)
    {
      for (auto & i : p)
	{
	  parents.push_back({i.first.clone(),i.second});
	}
    }
    
    virtual ~vnode()
    {}
    
  public:
    virtual std::shared_ptr<node<T_p,T_var1> > get_node()const
    {
      throw mcmc_exception("should never be called");
      return std::shared_ptr<node<T_p,T_var1> >();
    }

  public:
    virtual std::shared_ptr<vnode<T_p,T_var1> > clone()const
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new vnode<T_p,T_var1>(*this));
    }

  public:
    eq_vnode<T_p,T_var1> operator()(size_t n)const
    {
      return eq_vnode<T_p,T_var1>({*this,n});
    }
  };

  template <typename T_p,typename T_var1>
  class eq_vnode
    :public vnode<T_p,T_var1>
  {
  public:
    eq_vnode(const std::pair<const vnode<T_p,T_var1>&,size_t>& p)
      :vnode<T_p,T_var1>("eq",std::string("eq")+node_count<eq_vnode<T_p,T_var1> >(),{p})
    {
      this->binded=true;
    }

  public:
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new eq_node<T_p,T_var1>);
    }

  public:
    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new eq_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  class eq_node
    :public deterministic_node<T_p,T_var1>
  {
  public:
    eq_node()
      :deterministic_node<T_p,T_var1>(1,1)
      {}

    T_var1 do_value(size_t idx,size_t obsid)const override
    {
      return this->parent(0,obsid);
    }
  };
  
  
  template <typename T_p,typename T_var1>
  class new_node
  {
  public:
    std::string name;
  public:
    new_node(const std::string& n)
      :name(n)
    {}
    
    auto operator=(const vnode<T_p,T_var1>& rhs)
    {
      const_cast<vnode<T_p,T_var1>&>(rhs).name=name;
      cout<<typeid(rhs).name() << " " << rhs.name<<endl;
      return rhs;
    }
  };
}

#endif
