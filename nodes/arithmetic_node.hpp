#ifndef ARITHMETIC_NODE
#define ARITHMETIC_NODE
#include <core/tp_aware_dtm_node.hpp>
#include <helper/node_counter.hpp>
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
}

#endif
