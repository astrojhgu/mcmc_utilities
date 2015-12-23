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
  template <typename T>
  class add_node
    :public deterministic_node<T>
  {
  public:
    add_node()
      :deterministic_node<T>(2,1)
    {}
    
    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]+parent[1];
    }
  };

  template <typename T>
  class add_node_factory
    :public abstract_node_factory<T>
  {
  public:
    add_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new add_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T>
  class add_vnode
    :public vnode<T>
  {
  public:
    add_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("add",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new add_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new add_vnode<T>(*this));
    }
  };
  
  template <typename T>
  add_vnode<T> operator+(const vnode<T>& n1,
		 const vnode<T>& n2)
  {
    auto result= add_vnode<T>(std::string("add")+node_count<add_vnode<T> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }
  
  /////sub////
  template <typename T>
  class sub_node
    :public deterministic_node<T>
  {
  public:
    sub_node()
      :deterministic_node<T>(2,1)
    {}
    
    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]-parent[1];
    }
  };

  template <typename T>
  class sub_node_factory
    :public abstract_node_factory<T>
  {
  public:
    sub_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new sub_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T>
  class sub_vnode
    :public vnode<T>
  {
  public:
    sub_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("sub",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new sub_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new sub_vnode<T>(*this));
    }
  };
  
  template <typename T>
  sub_vnode<T> operator-(const vnode<T>& n1,
		 const vnode<T>& n2)
  {
    auto result = sub_vnode<T>(std::string("sub")+node_count<sub_vnode<T> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }

  /////neg////
  template <typename T>
  class neg_node
    :public deterministic_node<T>
  {
  public:
    neg_node()
      :deterministic_node<T>(1,1)
    {}
    
    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return -parent[0];
    }
  };

  template <typename T>
  class neg_node_factory
    :public abstract_node_factory<T>
  {
  public:
    neg_node_factory()
      :abstract_node_factory<T>({"op1"},{"result"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new neg_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T>
  class neg_vnode
    :public vnode<T>
  {
  public:
    neg_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("neg",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new neg_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new neg_vnode<T>(*this));
    }
  };

    /////pos////
  template <typename T>
  class pos_node
    :public deterministic_node<T>
  {
  public:
    pos_node()
      :deterministic_node<T>(1,1)
    {}
    
    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0];
    }
  };

  template <typename T>
  class pos_node_factory
    :public abstract_node_factory<T>
  {
  public:
    pos_node_factory()
      :abstract_node_factory<T>({"op1"},{"result"},{})
    {}
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new pos_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T>
  class pos_vnode
    :public vnode<T>
  {
  public:
    pos_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("pos",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new pos_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new pos_vnode<T>(*this));
    }
  };


  /////mul////
  template <typename T>
  class mul_node
    :public deterministic_node<T>
  {
  public:
    mul_node()
      :deterministic_node<T>(2,1)
    {}
    
    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]*parent[1];
    }
  };
  
  template <typename T>
  class mul_node_factory
    :public abstract_node_factory<T>
  {
  public:
    mul_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new mul_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  template <typename T>
  class mul_vnode
    :public vnode<T>
  {
  public:
    mul_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("mul",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new mul_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new mul_vnode<T>(*this));
    }
  };
  
  template <typename T>
  mul_vnode<T> operator*(const vnode<T>& n1,
		 const vnode<T>& n2)
  {
    auto result= mul_vnode<T>(std::string("mul")+node_count<mul_vnode<T> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }
  
  /////div////
  template <typename T>
  class div_node
    :public deterministic_node<T>
  {
  public:
    div_node()
      :deterministic_node<T>(2,1)
    {}
    
    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]/parent[1];
    }
  };

  template <typename T>
  class div_node_factory
    :public abstract_node_factory<T>
  {
  public:
    div_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new div_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T>
  class div_vnode
    :public vnode<T>
  {
  public:
    div_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("div",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new div_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new div_vnode<T>(*this));
    }
  };
  
  template <typename T>
  div_vnode<T> operator/(const vnode<T>& n1,
		 const vnode<T>& n2)
  {
    auto result= div_vnode<T>(std::string("div")+node_count<div_vnode<T> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }
  
  /////pow////
  template <typename T>
  class pow_node
    :public deterministic_node<T>
  {
  public:
    pow_node()
      :deterministic_node<T>(2,1)
    {}
    
    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return std::pow(parent[0],parent[1]);
    }
  };

  template <typename T>
  class pow_node_factory
    :public abstract_node_factory<T>
  {
  public:
    pow_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new pow_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T>
  class pow_vnode
    :public vnode<T>
  {
  public:
    pow_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("pow",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new pow_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new pow_vnode<T>(*this));
    }
  };
  
  template <typename T>
  pow_vnode<T> operator/(const vnode<T>& n1,
		 const vnode<T>& n2)
  {
    auto result= pow_vnode<T>(std::string("pow")+node_count<pow_vnode<T> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }

  
}

#endif
