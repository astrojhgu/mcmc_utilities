#ifndef FUNC_NODE_HPP
#define FUNC_NODE_HPP

#include <core/deterministic_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T>
  class func_node
    :public deterministic_node<T>
  {
    std::function<T (const std::vector<T>&)> func;
    int n;
  public:
    //template <typename Tf>
    func_node(T (*f)(const std::vector<T>&),int n1)
      :deterministic_node<T>(n1,1),func(f),n(n1)
    {}

    func_node(const func_node<T>& rhs)
      :deterministic_node<T>(rhs.n,1),func(rhs.func),n(rhs.n)
    {}
    

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      return func(parent);
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      return std::shared_ptr<node<T> >(new func_node(*this));
    }
  };
}


#endif
