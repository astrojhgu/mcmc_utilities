#ifndef ARITHMETIC_NODE
#define ARITHMETIC_NODE
#include <core/tp_aware_dtm_node.hpp>
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
  template <typename T,template <typename TE> class T_vector>
  class add_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    add_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]+parent[1];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new add_node);
    }
    
    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      order o1=this->get_parent_order(0,pn,n);
      order o2=this->get_parent_order(1,pn,n);
      if(o1.n>1||o1.n<0||
	 !o1.poly||
	 o2.n>1||o1.n<0||
	 !o2.poly)
	{
	  return order{0,false,false};
	}
      return order{std::max(o1.n,o2.n),(o1.n==o2.n&&o1.homo&&o2.homo),true};
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class add_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    add_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new add_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T,template <typename TE> class T_vector>
  class add_vnode
    :public vnode<T,T_vector>
  {
  public:
    add_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T,T_vector>&,size_t> >& p)
      :vnode<T,T_vector>("add",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new add_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new add_vnode<T,T_vector>(*this));
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  add_vnode<T,T_vector> operator+(const vnode<T,T_vector>& n1,
		 const vnode<T,T_vector>& n2)
  {
    auto result= add_vnode<T,T_vector>(std::string("add")+node_count<add_vnode<T,T_vector> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }
  
  /////sub////
  template <typename T,template <typename TE> class T_vector>
  class sub_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    sub_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]-parent[1];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new sub_node);
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      order o1=this->get_parent_order(0,pn,n);
      order o2=this->get_parent_order(1,pn,n);
      if(o1.n>1||o1.n<0||
	 !o1.poly||
	 o2.n>1||o1.n<0||
	 !o2.poly)
	{
	  return order{0,false,false};
	}
      return order{std::max(o1.n,o2.n),(o1.n==o2.n&&o1.homo&&o2.homo),true};
    }

  };

  template <typename T,template <typename TE> class T_vector>
  class sub_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    sub_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new sub_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T,template <typename TE> class T_vector>
  class sub_vnode
    :public vnode<T,T_vector>
  {
  public:
    sub_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T,T_vector>&,size_t> >& p)
      :vnode<T,T_vector>("sub",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new sub_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new sub_vnode<T,T_vector>(*this));
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  sub_vnode<T,T_vector> operator-(const vnode<T,T_vector>& n1,
		 const vnode<T,T_vector>& n2)
  {
    auto result = sub_vnode<T,T_vector>(std::string("sub")+node_count<sub_vnode<T,T_vector> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }

  /////neg////
  template <typename T,template <typename TE> class T_vector>
  class neg_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    neg_node()
      :tp_aware_dtm_node<T,T_vector>(1,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return -parent[0];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new neg_node);
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      order o=this->get_parent_order(0,pn,n);
      return o;
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class neg_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    neg_node_factory()
      :abstract_node_factory<T,T_vector>({"op1"},{"result"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new neg_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T,template <typename TE> class T_vector>
  class neg_vnode
    :public vnode<T,T_vector>
  {
  public:
    neg_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T,T_vector>&,size_t> >& p)
      :vnode<T,T_vector>("neg",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new neg_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new neg_vnode<T,T_vector>(*this));
    }
  };

    /////pos////
  template <typename T,template <typename TE> class T_vector>
  class pos_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    pos_node()
      :tp_aware_dtm_node<T,T_vector>(1,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new pos_node);
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      order o=this->get_parent_order(0,pn,n);
      return o;
    }

  };

  template <typename T,template <typename TE> class T_vector>
  class pos_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    pos_node_factory()
      :abstract_node_factory<T,T_vector>({"op1"},{"result"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new pos_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T,template <typename TE> class T_vector>
  class pos_vnode
    :public vnode<T,T_vector>
  {
  public:
    pos_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T,T_vector>&,size_t> >& p)
      :vnode<T,T_vector>("pos",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new pos_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new pos_vnode<T,T_vector>(*this));
    }
  };


  /////mul////
  template <typename T,template <typename TE> class T_vector>
  class mul_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    mul_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]*parent[1];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new mul_node);
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      order o1=this->get_parent_order(0,pn,n);
      order o2=this->get_parent_order(1,pn,n);
      
      if(!o1.poly||!o2.poly)
	{
	  return order{0,false,false};
	}

      return order{o1.n+o2.n,(o1.homo&&o2.homo),true};
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  class mul_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    mul_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new mul_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  template <typename T,template <typename TE> class T_vector>
  class mul_vnode
    :public vnode<T,T_vector>
  {
  public:
    mul_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T,T_vector>&,size_t> >& p)
      :vnode<T,T_vector>("mul",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new mul_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new mul_vnode<T,T_vector>(*this));
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  mul_vnode<T,T_vector> operator*(const vnode<T,T_vector>& n1,
		 const vnode<T,T_vector>& n2)
  {
    auto result= mul_vnode<T,T_vector>(std::string("mul")+node_count<mul_vnode<T,T_vector> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }
  
  /////div////
  template <typename T,template <typename TE> class T_vector>
  class div_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    div_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]/parent[1];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new div_node);
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      order o1=this->get_parent_order(0,pn,n);
      order o2=this->get_parent_order(1,pn,n);
      
      if(!o1.poly||!o2.poly)
	{
	  return order{0,false,false};
	}

      if(o1.homo&&o2.homo)
	{
	  return order{o1.n-o2.n,true,true};
	}
      else
	{
	  return order{o1.n-o2.n,false,false};
	}
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class div_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    div_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new div_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
  
  template <typename T,template <typename TE> class T_vector>
  class div_vnode
    :public vnode<T,T_vector>
  {
  public:
    div_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T,T_vector>&,size_t> >& p)
      :vnode<T,T_vector>("div",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new div_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new div_vnode<T,T_vector>(*this));
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  div_vnode<T,T_vector> operator/(const vnode<T,T_vector>& n1,
		 const vnode<T,T_vector>& n2)
  {
    auto result= div_vnode<T,T_vector>(std::string("div")+node_count<div_vnode<T,T_vector> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }
  
  /////pow////
  template <typename T,template <typename TE> class T_vector>
  class pow_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    pow_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return std::pow(parent[0],parent[1]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new pow_node);
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      order o1=this->get_parent_order(0,pn,n);
      order o2=this->get_parent_order(1,pn,n);
      
      if(!o1.poly||o1.n!=0||
	 !o2.poly||o2.n!=0)
	{
	  return order{0,false,false};
	}
      return order{0,true,true};
    }
    
  };

  template <typename T,template <typename TE> class T_vector>
  class pow_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    pow_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new pow_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };
  
  
  template <typename T,template <typename TE> class T_vector>
  class pow_vnode
    :public vnode<T,T_vector>
  {
  public:
    pow_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T,T_vector>&,size_t> >& p)
      :vnode<T,T_vector>("pow",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new pow_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new pow_vnode<T,T_vector>(*this));
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  pow_vnode<T,T_vector> operator/(const vnode<T,T_vector>& n1,
		 const vnode<T,T_vector>& n2)
  {
    auto result= pow_vnode<T,T_vector>(std::string("pow")+node_count<pow_vnode<T,T_vector> >(),{{n1,0},{n2,0}});
    result.named=false;
    return result;
  }

  
}

#endif
