#ifndef POISSON_NODE_HPP
#define POISSON_NODE_HPP
#include <helper/node_counter.hpp>
#include <math/distributions.hpp>
#include <math/functions.hpp>
#include <core/forward_sampleable_node.hpp>
#include <helper/abstract_node_factory.hpp>
namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class poisson_node
    :public forward_sampleable_node<T,T_vector>
  {
  public:
    poisson_node()
      :forward_sampleable_node<T,T_vector>(1,1)
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

    T do_regulate(size_t idx,const T& x)const override
    {
      return std::floor(x);
    }

    T_vector<T> do_candidate_points()const override
    {
      /*
      T_vector<T> result((int)(this->parent(0))*10+1);
      for(int i=0;i<result.size();++i)
	{
	  result[i]=i;
	}
      return result;
      */
      return T_vector<T>();
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,this->parent(0));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new poisson_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };

  template <typename T,template <typename TE> class T_vector>
  class poisson_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    poisson_node_factory()
      :abstract_node_factory<T,T_vector>({"lambda"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new poisson_node<T,T_vector>());
    }


    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };
}


#endif
