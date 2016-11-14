#ifndef UNARY_NODE_HPP
#define UNARY_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <functional>


namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class unary_node
    :public cached_dtm_node<T,T_vector>
  {
  private:
    std::function<T (const T&)> func;
  public:
    unary_node(const std::function<T (const T&)>& f)
      :cached_dtm_node<T,T_vector>(1,1),func(f)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return func(parent[0]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new unary_node(func));
    }

  };
}


#endif
