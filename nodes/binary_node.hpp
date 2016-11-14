#ifndef BINARY_NODE_HPP
#define BINARY_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>
#include <functional>


namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class binary_node
    :public cached_dtm_node<T,T_vector>
  {
  private:
    std::function<T (const T&)> func;
  public:
    binary_node(const std::function<T (const T&,const T&)>& f)
      :cached_dtm_node<T,T_vector>(1,1),func(f)
    {}

    T do_calc(size_t idx,const std::vector<T,T_vector>& parent)const override
    {
      return func(this->parent[0]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new binary_node(alpha,beta);
      return std::shared_ptr<node<T,T_vector> >(p);
    }
    
  };
}


#endif
