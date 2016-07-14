#ifndef LT_NODE_HPP
#define LT_NODE_HPP

#include <core/tp_aware_dtm_node.hpp>
#include <core/error_handler.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class lt_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    lt_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]<parent[1]?1:0;
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new lt_node();
      return std::shared_ptr<node<T,T_vector> >(p);
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
  class lt_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    lt_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new lt_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class le_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    le_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]<=parent[1]?1:0;
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new le_node();
      return std::shared_ptr<node<T,T_vector> >(p);
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
  class le_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    le_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new le_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };


  template <typename T,template <typename TE> class T_vector>
  class gt_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    gt_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]>parent[1]?1:0;
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new gt_node();
      return std::shared_ptr<node<T,T_vector> >(p);
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
  class gt_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    gt_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new gt_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class ge_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    ge_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]>=parent[1]?1:0;
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new ge_node();
      return std::shared_ptr<node<T,T_vector> >(p);
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
  class ge_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    ge_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new ge_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class eq_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    eq_node()
      :tp_aware_dtm_node<T,T_vector>(2,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]==parent[1]?1:0;
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new eq_node();
      return std::shared_ptr<node<T,T_vector> >(p);
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
  class eq_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    eq_node_factory()
      :abstract_node_factory<T,T_vector>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new eq_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };

}


#endif
