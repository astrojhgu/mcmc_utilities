#ifndef FUNC_NODE_HPP
#define FUNC_NODE_HPP

#include <core/cached_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class func_node
    :public cached_dtm_node<T,T_vector>
  {
    std::function<T (const T_vector<T>&)> func;
    int n;
  public:
    //template <typename Tf>
    func_node(T (*f)(const T_vector<T>&),int n1)
      :cached_dtm_node<T,T_vector>(n1,1),func(f),n(n1)
    {}

    func_node(const func_node<T,T_vector>& rhs)
      :cached_dtm_node<T,T_vector>(rhs.n,1),func(rhs.func),n(rhs.n)
    {}
    

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return func(parent);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new func_node(*this));
    }
  };
}


#endif
