#ifndef CHILDREN_SELECTOR_NODE_HPP
#define CHILDREN_SELECTOR_NODE_HPP

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
  class children_selector_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    children_selector_node(int noutput)
      :tp_aware_dtm_node<T,T_vector>(2,noutput)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return parent.at(1);
    }

    std::set<stochastic_node<T,T_vector>* > do_enumerate_stochastic_children()const override
    {
      std::set<stochastic_node<T,T_vector>* > result;
      size_t n=static_cast<size_t>(this->parent(0));
      std::for_each(std::begin(this->stochastic_children),std::end(this->stochastic_children),[&result,n](const std::pair<stochastic_node<T,T_vector>*,size_t>& i)
		   {
		     if(i.second==n)
		       {
			 result.insert(i.first);
		       }
		   });
      for(auto& p:this->deterministic_children)
	{
	  if(p.second!=n)
	    {
	      continue;
	    }
	  auto c=p.first->enumerate_stochastic_children();
	  result.insert(std::begin(c),std::end(c));
	}
      return std::move(result);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new children_selector_node(this->num_of_dims());
      return std::shared_ptr<node<T,T_vector> >(p);
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      for(size_t i=0;i<this->num_of_parents();++i)
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
