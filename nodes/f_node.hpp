#ifndef F_NODE_HPP
#define F_NODE_HPP
#include <core/forward_sampleable_node.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class f_node
    :public forward_sampleable_node<T,T_vector>
  {
  private:
  public:
    f_node()
      :forward_sampleable_node<T,T_vector>(2,std::numeric_limits<T>::epsilon())
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      T n=this->parent(0);
      T m=this->parent(1);
      return lgamma((m+n)/2)-lgamma(n/2)-lgamma(m/2)+(n/2)*std::log(n/m)+(n/2-1)*std::log(x)-(m+n)/2*(1+n*x/m);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      return std::make_pair(std::numeric_limits<T>::epsilon(),100);
    }

    void do_init_value(size_t) override
    {
      T n=this->parent(0);
      T m=this->parent(1);
      this->set_value(0,m);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new f_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };
  
  template <typename T,template <typename TE> class T_vector>
  class f_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    f_node_factory()
      :abstract_node_factory<T,T_vector>({"n","m"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new f_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
}

#endif
