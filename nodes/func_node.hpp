#ifndef FUNC_NODE_HPP
#define FUNC_NODE_HPP

#include <core/tp_aware_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class func_node
    :public tp_aware_dtm_node<T,T_vector>
  {
    std::function<T (const T_vector<T>&)> func;
    int n;
  public:
    //template <typename Tf>
    func_node(T (*f)(const T_vector<T>&),int n1)
      :tp_aware_dtm_node<T,T_vector>(n1,1),func(f),n(n1)
    {}

    func_node(const func_node<T,T_vector>& rhs)
      :tp_aware_dtm_node<T,T_vector>(rhs.n,1),func(rhs.func),n(rhs.n)
    {}
    

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return func(parent);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new func_node(*this));
    }

    order do_get_order(const node<T,T_vector>* pn,int n)const override
    {
      for (int i=0;i<this->num_of_parents();++i)
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
}


#endif
