#ifndef SIN_NODE_HPP
#define SIN_NODE_HPP

#include <core/differentiable_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class sin_node
    :public differentiable_dtm_node<T,T_vector>
  {
  public:
    sin_node()
      :differentiable_dtm_node<T,T_vector>(1,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return std::sin(parent[0]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new sin_node;
      return std::shared_ptr<node<T,T_vector> >(p);
    }
    
  };


  template <typename T,template <typename TE> class T_vector>
  class sin_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    sin_node_factory()
      :abstract_node_factory<T,T_vector>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new sin_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };
}


#endif
