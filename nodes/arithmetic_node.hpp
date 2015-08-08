#ifndef ARITHMETIC_NODE
#define ARITHMETIC_NODE
#include <core/deterministic_node.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <initializer_list>

namespace mcmc_utilities
{
  /////add////
  template <typename T_p,typename T_var1>
  class add_node
    :public deterministic_node<T_p,T_var1>
  {
  private:
    
  public:
    add_node()
      :deterministic_node<T_p,T_var1>(2,1)
    {}
    
    T_var1 do_value(size_t idx,size_t obsid)const override
    {
      return this->parent(0,obsid)+this->parent(1,obsid);
    }
  };
  
  
  
  template <typename T_p,typename T_var1>
  class _add_vnode
    :public _vnode<T_p, T_var1>
  {
  public:
    _add_vnode(std::string n,const std::initializer_list<std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> >& p)
      :_vnode<T_p,T_var1>("add",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new add_node<T_p,T_var1>);
    }
  };
  
  template <typename T_p,typename T_var1>
  auto operator+(const cpnt<T_p,T_var1>& n1,
		 const cpnt<T_p,T_var1>& n2)
  {
    return cpnt<T_p,T_var1>(std::shared_ptr<_vnode<T_p,T_var1> >(
								 new _add_vnode<T_p,T_var1>(std::string("add")+node_count<_add_vnode<T_p,T_var1> >(),{{n1.pn,n1.n},{n2.pn,n2.n}})
								 ));
  }


  /////sub////
  template <typename T_p,typename T_var1>
  class sub_node
    :public deterministic_node<T_p,T_var1>
  {
  private:
    
  public:
    sub_node()
      :deterministic_node<T_p,T_var1>(2,1)
    {}
    
    T_var1 do_value(size_t idx,size_t obsid)const override
    {
      return this->parent(0,obsid)-this->parent(1,obsid);
    }
  };
  
  
  
  template <typename T_p,typename T_var1>
  class _sub_vnode
    :public _vnode<T_p, T_var1>
  {
  public:
    _sub_vnode(std::string n,const std::initializer_list<std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> >& p)
      :_vnode<T_p,T_var1>("sub",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new sub_node<T_p,T_var1>);
    }
  };
  
  template <typename T_p,typename T_var1>
  auto operator-(const cpnt<T_p,T_var1>& n1,
		 const cpnt<T_p,T_var1>& n2)
  {
    return cpnt<T_p,T_var1>(std::shared_ptr<_vnode<T_p,T_var1> >(
								 new _sub_vnode<T_p,T_var1>(std::string("sub")+node_count<_sub_vnode<T_p,T_var1> >(),{{n1.pn,n1.n},{n2.pn,n2.n}})
								 ));
  }


    /////mul////
  template <typename T_p,typename T_var1>
  class mul_node
    :public deterministic_node<T_p,T_var1>
  {
  private:
    
  public:
    mul_node()
      :deterministic_node<T_p,T_var1>(2,1)
    {}
    
    T_var1 do_value(size_t idx,size_t obsid)const override
    {
      return this->parent(0,obsid)*this->parent(1,obsid);
    }
  };
  
  
  
  template <typename T_p,typename T_var1>
  class _mul_vnode
    :public _vnode<T_p, T_var1>
  {
  public:
    _mul_vnode(std::string n,const std::initializer_list<std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> >& p)
      :_vnode<T_p,T_var1>("mul",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new mul_node<T_p,T_var1>);
    }
  };
  
  template <typename T_p,typename T_var1>
  auto operator*(const cpnt<T_p,T_var1>& n1,
		 const cpnt<T_p,T_var1>& n2)
  {
    return cpnt<T_p,T_var1>(std::shared_ptr<_vnode<T_p,T_var1> >(
								 new _mul_vnode<T_p,T_var1>(std::string("mul")+node_count<_mul_vnode<T_p,T_var1> >(),{{n1.pn,n1.n},{n2.pn,n2.n}})
								 ));
  }

      /////div////
  template <typename T_p,typename T_var1>
  class div_node
    :public deterministic_node<T_p,T_var1>
  {
  private:
    
  public:
    div_node()
      :deterministic_node<T_p,T_var1>(2,1)
    {}
    
    T_var1 do_value(size_t idx,size_t obsid)const override
    {
      return this->parent(0,obsid)/this->parent(1,obsid);
    }
  };
  
  
  
  template <typename T_p,typename T_var1>
  class _div_vnode
    :public _vnode<T_p, T_var1>
  {
  public:
    _div_vnode(std::string n,const std::initializer_list<std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> >& p)
      :_vnode<T_p,T_var1>("div",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new div_node<T_p,T_var1>);
    }
  };
  
  template <typename T_p,typename T_var1>
  auto operator/(const cpnt<T_p,T_var1>& n1,
		 const cpnt<T_p,T_var1>& n2)
  {
    return cpnt<T_p,T_var1>(std::shared_ptr<_vnode<T_p,T_var1> >(
								 new _div_vnode<T_p,T_var1>(std::string("div")+node_count<_div_vnode<T_p,T_var1> >(),{{n1.pn,n1.n},{n2.pn,n2.n}})
								 ));
  }
}

#endif
