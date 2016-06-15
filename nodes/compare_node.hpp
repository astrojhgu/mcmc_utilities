#ifndef LT_NODE_HPP
#define LT_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <core/error_handler.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class lt_node
    :public cached_dtm_node<T>
  {
  public:
    lt_node()
      :cached_dtm_node<T>(2,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]<parent[1]?1:0;
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new lt_node();
      return std::shared_ptr<node<T> >(p);
    }
  };

  template <typename T>
  class lt_node_factory
    :public abstract_node_factory<T>
  {
  public:
    lt_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new lt_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };

  template <typename T>
  class le_node
    :public cached_dtm_node<T>
  {
  public:
    le_node()
      :cached_dtm_node<T>(2,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]<=parent[1]?1:0;
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new le_node();
      return std::shared_ptr<node<T> >(p);
    }
  };

  template <typename T>
  class le_node_factory
    :public abstract_node_factory<T>
  {
  public:
    le_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new le_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };


  template <typename T>
  class gt_node
    :public cached_dtm_node<T>
  {
  public:
    gt_node()
      :cached_dtm_node<T>(2,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]>parent[1]?1:0;
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new gt_node();
      return std::shared_ptr<node<T> >(p);
    }
  };


  template <typename T>
  class gt_node_factory
    :public abstract_node_factory<T>
  {
  public:
    gt_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new gt_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };

  template <typename T>
  class ge_node
    :public cached_dtm_node<T>
  {
  public:
    ge_node()
      :cached_dtm_node<T>(2,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]>=parent[1]?1:0;
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new ge_node();
      return std::shared_ptr<node<T> >(p);
    }
  };

  template <typename T>
  class ge_node_factory
    :public abstract_node_factory<T>
  {
  public:
    ge_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new ge_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };

  template <typename T>
  class eq_node
    :public cached_dtm_node<T>
  {
  public:
    eq_node()
      :cached_dtm_node<T>(2,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]==parent[1]?1:0;
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new eq_node();
      return std::shared_ptr<node<T> >(p);
    }
  };

  template <typename T>
  class eq_node_factory
    :public abstract_node_factory<T>
  {
  public:
    eq_node_factory()
      :abstract_node_factory<T>({"op1","op2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new eq_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }
  };

}


#endif
