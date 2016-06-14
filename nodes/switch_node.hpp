#ifndef SWITCH_NODE_HPP
#define SWITCH_NODE_HPP

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
  class switch_node
    :public cached_dtm_node<T>
  {
  public:
    switch_node(int ninputs)
      :cached_dtm_node<T>(ninputs+1,1)
    {}

    T do_calc(size_t idx,const std::vector<T>& parent)const override
    {
      int n=size_t(parent.back());
      if(n<0||n>=this->num_of_parents()-1)
	{
	  throw var_out_of_range();
	}
      return parent.at(n);
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new switch_node(this->num_of_parents()-1);
      return std::shared_ptr<node<T> >(p);
    }
  };
}


#endif
