#ifndef BIN_NODE_HPP
#define BIN_NODE_HPP
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <math/distributions.hpp>
#include <math/functions.hpp>
#include <core/stochastic_node.hpp>
namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class bin_node
    :public stochastic_node<T_p,T_var1>
  {
  public:
    bin_node()
      :stochastic_node<T_p,T_var1>(2,1)
    {}

  public:
    T_p do_log_prob()const override
    {
      T_p result=logdbin(this->value(0),this->parent(0),this->parent(1));
      return result;
    }

    std::pair<T_var1,T_var1> do_var_range()const override
    {
      return std::pair<T_var1,T_var1>(0,this->parent(1));
    }

    bool is_continuous(size_t idx)const override
    {
      if(idx==0)
	{
	  return false;
	}
      else
	{
	  return true;
	}
    }

    std::vector<T_var1> do_candidate_points()const override
    {
      std::vector<T_var1> result((int)(this->parent(1))+1);
      for(int i=0;i<result.size();++i)
	{
	  result[i]=i;
	}
      return result;
    }
  };

  template <typename T_p,typename T_var1>
  class bin_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    bin_node_factory()
      :abstract_node_factory<T_p,T_var1>({"p","n"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(const std::vector<T_var1>& hparam)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new bin_node<T_p,T_var1>());
    }


    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };

  template <typename T_p,typename T_var1>
  class bin_vnode
    :public vnode<T_p,T_var1>
  {
  private:
    
  public:
    
    bin_vnode(std::string n,
		  const std::pair<const vnode<T_p,T_var1>&,size_t>& p1,
		  const std::pair<const vnode<T_p,T_var1>&,size_t>& p2)
      :vnode<T_p,T_var1>("bin",n,{p1,p2})
    {
      this->binded=true;
    }

    bin_vnode(const std::pair<const vnode<T_p,T_var1>&,size_t>& p1,
	      const std::pair<const vnode<T_p,T_var1>&,size_t>& p2)
      :vnode<T_p,T_var1>("bin",std::string("bin")+node_count<bin_vnode<T_p,T_var1> >(),{p1,p2})
    {
      this->binded=true;
      this->named=false;
    }
    
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new bin_node<T_p,T_var1>());
    }
    
    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new bin_vnode<T_p,T_var1>(*this));
    }
  };
  
  using vbin=bin_vnode<double,double>;
}


#endif
