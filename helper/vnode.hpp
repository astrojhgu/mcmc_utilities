#ifndef VNODE_HPP
#define VNODE_HPP
#include <string>
#include <vector>
#include <initializer_list>
#include <map>
#include <memory>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class _vnode
  {
  public:
    std::string type;
    std::string name;
    std::vector<std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> > parents;
    bool binded;

    bool added;
    
    _vnode(const std::string& n)
      :_vnode("unbinded",n,{})
    {}
    
    _vnode(const std::string& t,const std::string& n)
      :_vnode(t,n,{})
    {}
    
    _vnode(const std::string& t,const std::string& n,const std::initializer_list<std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> >& p)
      :type(t),name(n),parents(p),binded(false),added(false)
    {
    }
    
    virtual ~_vnode()
    {}
    
  public:
    virtual std::shared_ptr<node<T_p,T_var1> > get_node()const
    {
      throw mcmc_exception("should never be called");
      return std::shared_ptr<node<T_p,T_var1> >();
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
    
    auto operator=(const std::shared_ptr<_vnode<T_p,T_var1> >& rhs)
    {
      rhs->name=name;
      return rhs;
    }
  };
  
  
  template <typename T_p,typename T_var1>
  auto vnode(std::string n)
  {
    return std::shared_ptr<_vnode<T_p,T_var1> >(new _vnode<T_p,T_var1>(n));
  }
}

#endif
