#ifndef COS_NODE_HPP
#define COS_NODE_HPP

#include <core/tp_aware_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class cos_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    cos_node()
      :tp_aware_dtm_node<T,T_vector>(1,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return std::cos(parent[0]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new cos_node);
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      for (size_t i=0;i<this->num_of_parents();++i)
	{
	  order o=this->get_parent_order(i,pn,n);
	  if(o.n!=0||
	     !o.poly)
	    {
	      return order{0,false,false};
	    }
	}
      return order{0,true,true};
    }

  };
  
  template <typename T,template <typename TE> class T_vector>
  class cos_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    cos_node_factory()
      :abstract_node_factory<T,T_vector>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new cos_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };
}


#endif
