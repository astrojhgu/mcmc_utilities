#ifndef BIN_NODE_HPP
#define BIN_NODE_HPP
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <math/distributions.hpp>
#include <math/functions.hpp>
#include <core/stochastic_node.hpp>
namespace mcmc_utilities
{
  template <typename T>
  class bin_node
    :public stochastic_node<T>
  {
  public:
    bin_node()
      :stochastic_node<T>(2,1)
    {}

  public:
    T do_log_prob()const override
    {
      T result=logdbin(this->value(0),this->parent(0),this->parent(1));
      return result;
    }

    std::pair<T,T> do_var_range()const override
    {
      return std::pair<T,T>(0,this->parent(1));
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
      std::vector<T> result((int)(this->parent(1))+1);
      for(int i=0;i<result.size();++i)
	{
	  result[i]=i;
	}
      return result;
    }

    void do_initialize(size_t n) override
    {
      this->set_value(0,(this->parent(0))*(this->parent(1)));
    }

    std::shared_ptr<node<T> > do_clone()const override
    {
      auto p=new bin_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T> >(p);
    }
  };

  template <typename T>
  class bin_node_factory
    :public abstract_node_factory<T>
  {
  public:
    bin_node_factory()
      :abstract_node_factory<T>({"p","n"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new bin_node<T>());
    }


    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };

  template <typename T>
  class bin_vnode
    :public vnode<T>
  {
  private:
    
  public:
    
    bin_vnode(std::string n,
		  const std::pair<const vnode<T>&,size_t>& p1,
		  const std::pair<const vnode<T>&,size_t>& p2)
      :vnode<T>("bin",n,{p1,p2})
    {
      this->binded=true;
    }

    bin_vnode(const std::pair<const vnode<T>&,size_t>& p1,
	      const std::pair<const vnode<T>&,size_t>& p2)
      :vnode<T>("bin",std::string("bin")+node_count<bin_vnode<T> >(),{p1,p2})
    {
      this->binded=true;
      this->named=false;
    }
    
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new bin_node<T>());
    }
    
    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new bin_vnode<T>(*this));
    }
  };
  
  using vbin=bin_vnode<double>;
}


#endif
