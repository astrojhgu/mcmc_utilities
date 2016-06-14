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
  template <typename T>
  class id_vnode;

  template <typename T>
  class id_node;
  
  template <typename T>
  class vnode
  {
  public:
    std::string type;
    std::string name;
    std::vector<std::pair<std::shared_ptr<vnode<T> >,size_t> > parents;
    bool binded;
    bool added;
    bool named;
    
    vnode(const std::string& n)
      :vnode("unbinded",n)
    {}
    
    vnode(const std::string& t,const std::string& n)
      :type(t),name(n),binded(false),added(false),named(true)
    {}
    
    vnode(const std::string& t,const std::string& n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :type(t),name(n),binded(false),added(false),named(true)
    {
      for (auto & i : p)
	{
	  parents.push_back({i.first.clone(),i.second});
	}
    }

    
    virtual ~vnode()
    {}
    
  public:
    virtual std::shared_ptr<node<T> > get_node()const
    {
      throw mcmc_exception("should never be called");
      return std::shared_ptr<node<T> >();
    }

  public:
    virtual std::shared_ptr<vnode<T> > clone()const
    {
      return std::shared_ptr<vnode<T> >(new vnode<T>(*this));
    }

    void set_name(const std::string& n)
    {
      name=n;
      named=true;
    }

  public:
    id_vnode<T> operator()(size_t n)const
    {
      return id_vnode<T>({*this,n});
    }
  };

  template <typename T>
  class id_vnode
    :public vnode<T>
  {
  public:
    id_vnode(const std::pair<const vnode<T>&,size_t>& p)
      :vnode<T>("eq",std::string("eq")+node_count<id_vnode<T> >(),{p})
    {
      this->binded=true;
      this->named=false;
    }

  public:
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new id_node<T>);
    }

  public:
    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new id_vnode<T>(*this));
    }
  };

  template <typename T>
  class id_node
    :public deterministic_node<T>
  {
  public:
    id_node()
      :deterministic_node<T>(1,1)
      {}

    T do_value(size_t idx)const override
    {
      return this->parent(0);
    }
  };
  
 
  using vn=vnode<double>;
}

#endif
