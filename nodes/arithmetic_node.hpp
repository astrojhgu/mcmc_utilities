#ifndef ARITHMETIC_NODE
#define ARITHMETIC_NODE
#include <core/deterministic_node.hpp>
#include <helper/node_counter.hpp>
#include <helper/vnode.hpp>
#include <memory>
#include <utility>
#include <string>
#include <initializer_list>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  /////add////
  template <typename T_p,typename T_var1>
  class add_node
    :public deterministic_node<T_p,T_var1>
  {
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
  class add_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    add_node_factory()
      :abstract_node_factory<T_p,T_var1>({"op1","op2"},{"result"},{},{})
    {}
    
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(
	     const std::vector<T_var1>& scalar_param,
	     const std::vector<std::vector<T_var1> >& vector_param)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new add_node<T_p,T_var1>);
    }      
  };
  
  
  template <typename T_p,typename T_var1>
  class add_vnode
    :public vnode<T_p, T_var1>
  {
  public:
    add_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T_p,T_var1>&,size_t> >& p)
      :vnode<T_p,T_var1>("add",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new add_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new add_vnode<T_p,T_var1>(*this));
    }
  };
  
  template <typename T_p,typename T_var1>
  add_vnode<T_p,T_var1> operator+(const vnode<T_p,T_var1>& n1,
		 const vnode<T_p,T_var1>& n2)
  {
    auto result= add_vnode<T_p,T_var1>(std::string("add")+node_count<add_vnode<T_p,T_var1> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }
  
  /////sub////
  template <typename T_p,typename T_var1>
  class sub_node
    :public deterministic_node<T_p,T_var1>
  {
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
  class sub_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    sub_node_factory()
      :abstract_node_factory<T_p,T_var1>({"op1","op2"},{"result"},{},{})
    {}
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(
	     const std::vector<T_var1>& scalar_param,
	     const std::vector<std::vector<T_var1> >& vector_param)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new sub_node<T_p,T_var1>);
    }      
  };
  
  
  template <typename T_p,typename T_var1>
  class sub_vnode
    :public vnode<T_p, T_var1>
  {
  public:
    sub_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T_p,T_var1>&,size_t> >& p)
      :vnode<T_p,T_var1>("sub",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new sub_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new sub_vnode<T_p,T_var1>(*this));
    }
  };
  
  template <typename T_p,typename T_var1>
  sub_vnode<T_p,T_var1> operator-(const vnode<T_p,T_var1>& n1,
		 const vnode<T_p,T_var1>& n2)
  {
    auto result = sub_vnode<T_p,T_var1>(std::string("sub")+node_count<sub_vnode<T_p,T_var1> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }


  /////mul////
  template <typename T_p,typename T_var1>
  class mul_node
    :public deterministic_node<T_p,T_var1>
  {
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
  class mul_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    mul_node_factory()
      :abstract_node_factory<T_p,T_var1>({"op1","op2"},{"result"},{},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(
		const std::vector<T_var1>& scalar_param,
		const std::vector<std::vector<T_var1> >& vector_param)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new mul_node<T_p,T_var1>);
    }      
  };
  
  template <typename T_p,typename T_var1>
  class mul_vnode
    :public vnode<T_p, T_var1>
  {
  public:
    mul_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T_p,T_var1>&,size_t> >& p)
      :vnode<T_p,T_var1>("mul",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new mul_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new mul_vnode<T_p,T_var1>(*this));
    }
  };
  
  template <typename T_p,typename T_var1>
  mul_vnode<T_p,T_var1> operator*(const vnode<T_p,T_var1>& n1,
		 const vnode<T_p,T_var1>& n2)
  {
    auto result= mul_vnode<T_p,T_var1>(std::string("mul")+node_count<mul_vnode<T_p,T_var1> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }
  
  /////div////
  template <typename T_p,typename T_var1>
  class div_node
    :public deterministic_node<T_p,T_var1>
  {
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
  class div_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    div_node_factory()
      :abstract_node_factory<T_p,T_var1>({"op1","op2"},{"result"},{},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(
	     const std::vector<T_var1>& scalar_param,
	     const std::vector<std::vector<T_var1> >& vector_param)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new div_node<T_p,T_var1>);
    }      
  };
  
  
  template <typename T_p,typename T_var1>
  class div_vnode
    :public vnode<T_p, T_var1>
  {
  public:
    div_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T_p,T_var1>&,size_t> >& p)
      :vnode<T_p,T_var1>("div",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new div_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new div_vnode<T_p,T_var1>(*this));
    }
  };
  
  template <typename T_p,typename T_var1>
  div_vnode<T_p,T_var1> operator/(const vnode<T_p,T_var1>& n1,
		 const vnode<T_p,T_var1>& n2)
  {
    auto result= div_vnode<T_p,T_var1>(std::string("div")+node_count<div_vnode<T_p,T_var1> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }
  
  
}

#endif
