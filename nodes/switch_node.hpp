#ifndef SWITCH_NODE_HPP
#define SWITCH_NODE_HPP

#include <core/tp_aware_dtm_node.hpp>
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
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    switch_node(int ninputs)
      :tp_aware_dtm_node<T,T_vector>(ninputs+1,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      int n=size_t(parent.back());
      if(n<0||n>=this->num_of_parents()-1)
	{
	  throw var_out_of_range();
	}
      return parent.at(n);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new switch_node(this->num_of_parents()-1);
      return std::shared_ptr<node<T,T_vector> >(p);
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      for(int i=0;i<this->num_of_parents();++i)
	{
	  order o=this->get_parent_order(0,pn,n);
	  //return order{0,false,false};
	  if(o.n!=0||
	     !o.poly)
	    {
	      return order{0,false,false};
	    }
	}
      return order{0,true,true};
    } 
  };
}


#endif
