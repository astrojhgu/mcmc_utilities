#ifndef COND_NODE_HPP
#define COND_NODE_HPP

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
  class cond_node
    :public cached_dtm_node<T>
  {
  public:
    cond_node()
      :cached_dtm_node<T>(3,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return parent[0]!=0?parent[1]:parent[2];
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new cond_node();
      return std::shared_ptr<node<T> >(p);
    }
  };

  template <typename T>
  class cond_node_factory
    :public abstract_node_factory<T>
  {
  public:
    cond_node_factory()
      :abstract_node_factory<T>({"s","in1","in2"},{"result"},{})
    {}
    
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new cond_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic");
    }

  };
  
}


#endif
