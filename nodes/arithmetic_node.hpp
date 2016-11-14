#ifndef ARITHMETIC_NODE
#define ARITHMETIC_NODE
#include <core/differentiable_dtm_node.hpp>
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
    :public differentiable_dtm_node<T,T_vector>
  {
  public:
    add_node()
      :differentiable_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]+parent[1];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new add_node);
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
    :public differentiable_dtm_node<T,T_vector>
  {
  public:
    sub_node()
      :differentiable_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]-parent[1];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new sub_node);
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
    :public differentiable_dtm_node<T,T_vector>
  {
  public:
    neg_node()
      :differentiable_dtm_node<T,T_vector>(1,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return -parent[0];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new neg_node);
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
    :public differentiable_dtm_node<T,T_vector>
  {
  public:
    pos_node()
      :differentiable_dtm_node<T,T_vector>(1,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new pos_node);
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
    :public differentiable_dtm_node<T,T_vector>
  {
  public:
    mul_node()
      :differentiable_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]*parent[1];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new mul_node);
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
    :public differentiable_dtm_node<T,T_vector>
  {
  public:
    div_node()
      :differentiable_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent[0]/parent[1];
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new div_node);
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
    :public differentiable_dtm_node<T,T_vector>
  {
  public:
    pow_node()
      :differentiable_dtm_node<T,T_vector>(2,1)
    {}
    
    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return std::pow(parent[0],parent[1]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new pow_node);
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
