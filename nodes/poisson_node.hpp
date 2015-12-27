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
  template <typename T>
  class poisson_node
    :public stochastic_node<T>
  {
  public:
    poisson_node()
      :stochastic_node<T>(1,1)
    {}

  public:
    T do_log_prob()const override
    {
      T result=logdpoisson(this->value(0),this->parent(0));
      return result;
    }

    std::pair<T,T> do_var_range()const override
    {
      return std::pair<T,T>(0,this->parent(0)*10);
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

    T do_regulate(const T& x)const override
    {
      return std::floor(x);
    }

    std::vector<T> do_candidate_points()const override
    {
      /*
      std::vector<T> result((int)(this->parent(0))*10+1);
      for(int i=0;i<result.size();++i)
	{
	  result[i]=i;
	}
      return result;
      */
      return std::vector<T>();
    }

    void do_initialize(size_t n)override
    {
      this->set_value(0,this->parent(0));
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new poisson_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T> >(p);
    }

  };

  template <typename T>
  class poisson_node_factory
    :public abstract_node_factory<T>
  {
  public:
    poisson_node_factory()
      :abstract_node_factory<T>({"lambda"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new poisson_node<T>());
    }


    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };

  template <typename T>
  class poisson_vnode
    :public vnode<T>
  {
  private:
    
  public:
    
    poisson_vnode(std::string n,
		  const std::pair<const vnode<T>&,size_t>& p1
		  )
      :vnode<T>("poisson",n,{p1})
    {
      this->binded=true;
    }

    poisson_vnode(const std::pair<const vnode<T>&,size_t>& p1
		  )
      :vnode<T>("poisson",std::string("poisson")+node_count<poisson_vnode<T> >(),{p1})
    {
      this->binded=true;
      this->named=false;
    }
    
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new poisson_node<T>());
    }
    
    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new poisson_vnode<T>(*this));
    }
  };
  
  using vpoisson=poisson_vnode<double>;
}


#endif
