#ifndef POISSON_NODE_HPP
#define POISSON_NODE_HPP
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <math/distributions.hpp>
#include <math/functions.hpp>
#include <core/stochastic_node.hpp>
#include <helper/abstract_node_factory.hpp>
namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class poisson_node
    :public stochastic_node<T_p,T_var1>
  {
  public:
    poisson_node()
      :stochastic_node<T_p,T_var1>(1,1)
    {}

  public:
    T_p do_log_prob()const override
    {
      T_p result=logdpoisson(this->value(0),this->parent(0));
      return result;
    }

    std::pair<T_var1,T_var1> do_var_range()const override
    {
      return std::pair<T_var1,T_var1>(0,this->parent(0)*10);
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
      std::vector<T_var1> result((int)(this->parent(0))*10+1);
      for(int i=0;i<result.size();++i)
	{
	  result[i]=i;
	}
      return result;
    }
  };

  template <typename T_p,typename T_var1>
  class poisson_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    poisson_node_factory()
      :abstract_node_factory<T_p,T_var1>({"lambda"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(const std::vector<T_var1>& hparam)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new poisson_node<T_p,T_var1>());
    }


    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };

  template <typename T_p,typename T_var1>
  class poisson_vnode
    :public vnode<T_p,T_var1>
  {
  private:
    
  public:
    
    poisson_vnode(std::string n,
		  const std::pair<const vnode<T_p,T_var1>&,size_t>& p1
		  )
      :vnode<T_p,T_var1>("poisson",n,{p1})
    {
      this->binded=true;
    }

    poisson_vnode(const std::pair<const vnode<T_p,T_var1>&,size_t>& p1
		  )
      :vnode<T_p,T_var1>("poisson",std::string("poisson")+node_count<poisson_vnode<T_p,T_var1> >(),{p1})
    {
      this->binded=true;
      this->named=false;
    }
    
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new poisson_node<T_p,T_var1>());
    }
    
    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new poisson_vnode<T_p,T_var1>(*this));
    }
  };
  
  using vpoisson=poisson_vnode<double,double>;
}


#endif
