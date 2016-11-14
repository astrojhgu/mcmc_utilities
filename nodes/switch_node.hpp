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
  template <typename T,template <typename TE> class T_vector>
  class switch_node
    :public cached_dtm_node<T,T_vector>
  {
  public:
    switch_node(int ninputs)
      :cached_dtm_node<T,T_vector>(ninputs+1,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      size_t n=size_t(parent.back());
      if(n+1>=this->num_of_parents())
	{
	  std::cerr<<"n="<<n<<std::endl;
	  std::cerr<<"nswitch="<<(this->num_of_parents()-1)<<std::endl;
	  throw var_out_of_range();
	}
      return parent.at(n);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new switch_node(this->num_of_parents()-1);
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };
}


#endif
