#ifndef LOG10_NODE_HPP
#define LOG10_NODE_HPP

#include <core/tp_aware_dtm_node.hpp>
#include <math/functions.hpp>
#include <helper/node_counter.hpp>
#include <memory>
#include <utility>
#include <string>



namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class log10_node
    :public tp_aware_dtm_node<T,T_vector>
  {
  public:
    log10_node()
      :tp_aware_dtm_node<T,T_vector>(1,1)
    {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return std::log10(parent[0]);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new log10_node);
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
  class log10_vnode
    :public vnode<T,T_vector>
  {
  public:
    log10_vnode(std::string n,
	      const std::pair<const vnode<T,T_vector>&,size_t>& p)
      :vnode<T,T_vector>("log10",n,{p})
    {
      this->binded=true;
    }

    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new log10_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new log10_vnode<T,T_vector>(*this));
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class log10_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    log10_node_factory()
      :abstract_node_factory<T,T_vector>({"x"},{"y"},{})
    {}
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new log10_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("deterministic node");
    }
  };

  template <typename T,template <typename TE> class T_vector>
  log10_vnode<T,T_vector> vlog10(const vnode<T,T_vector>& n1)
  {
    auto result= log10_vnode<T,T_vector>(std::string("log10")+node_count<log10_vnode<T,T_vector> >(),{n1,(size_t)0});
    result.named=false;
    return result;
  }

}


#endif
