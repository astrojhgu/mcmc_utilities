#ifndef BIN_NODE_HPP
#define BIN_NODE_HPP
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <math/distributions.hpp>
#include <math/functions.hpp>
#include <core/stochastic_node.hpp>
namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class bin_node
    :public stochastic_node<T,T_vector>
  {
  public:
    bin_node()
      :stochastic_node<T,T_vector>(2,1)
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

    T do_regulate(size_t idx,const T& x)const override
    {
      return std::floor(x);
    }

    T_vector<T> do_candidate_points()const override
    {
      T_vector<T> result((int)(this->parent(1))+1);
      for(unsigned int i=0;i<result.size();++i)
	{
	  result[i]=i;
	}
      return result;
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,(this->parent(0))*(this->parent(1)));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new bin_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class bin_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    bin_node_factory()
      :abstract_node_factory<T,T_vector>({"p","n"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new bin_node<T,T_vector>());
    }


    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class bin_vnode
    :public vnode<T,T_vector>
  {
  private:
    
  public:
    
    bin_vnode(std::string n,
		  const std::pair<const vnode<T,T_vector>&,size_t>& p1,
		  const std::pair<const vnode<T,T_vector>&,size_t>& p2)
      :vnode<T,T_vector>("bin",n,{p1,p2})
    {
      this->binded=true;
    }

    bin_vnode(const std::pair<const vnode<T,T_vector>&,size_t>& p1,
	      const std::pair<const vnode<T,T_vector>&,size_t>& p2)
      :vnode<T,T_vector>("bin",std::string("bin")+node_count<bin_vnode<T,T_vector> >(),{p1,p2})
    {
      this->binded=true;
      this->named=false;
    }
    
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new bin_node<T,T_vector>());
    }
    
    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new bin_vnode<T,T_vector>(*this));
    }
  };
  
  //using vbin=bin_vnode<double>;
}


#endif
