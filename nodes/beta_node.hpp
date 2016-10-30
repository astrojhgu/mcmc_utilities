#ifndef BETA_NODE_HPP
#define BETA_NODE_HPP
#include <core/forward_sampleable_node.hpp>
#include <helper/node_counter.hpp>
#include <limits>
#include <string>
#include <cmath>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class beta_node
    :public forward_sampleable_node<T,T_vector>
  {
  private:
    T alpha,beta;
  public:
    beta_node(T a,T b)
      :forward_sampleable_node<T,T_vector>(0,.5),alpha(a),beta(b)
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      return (alpha-1)*std::log(x)+(beta-1)*std::log(1-x)-lbeta(alpha,beta);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      return std::make_pair(std::numeric_limits<T>::epsilon(),T(1-std::numeric_limits<T>::epsilon()));
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,.5);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new beta_node(alpha,beta);
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  class beta_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    beta_node_factory()
      :abstract_node_factory<T,T_vector>({},{"x"},{"a","b"})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new beta_node<T,T_vector>(hparam[0],hparam[1]));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  //using vbeta=beta_vnode<double>;
}

#endif
