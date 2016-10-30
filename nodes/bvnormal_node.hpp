#ifndef BVNORMAL_NODE_HPP
#define BVNORMAL_NODE_HPP
#include <core/forward_sampleable_node.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class bvnormal_node
    :public forward_sampleable_node<T,T_vector>
  {
  private:
  public:
    bvnormal_node()
      :forward_sampleable_node<T,T_vector>(5,{0,0})
    {}
    
  private:
    T do_log_prob()const override
    {
      //const static T PI=std::atan(1.0)*4;
      T x1=this->value(0);
      T x2=this->value(1);
      T mu1=this->parent(0);
      T mu2=this->parent(1);
      T sigma1=this->parent(2);
      T sigma2=this->parent(3);
      T rho=this->parent(4);
      //return -(x-mu)*(x-mu)/(2*sigma*sigma)-std::log(sigma*std::sqrt(2*PI));
      T X1=(x1-mu1)/sigma1;
      T X2=(x2-mu2)/sigma2;
      T z=(X1*X1+X2*X2-2*rho*X1*X2)/(2*(1-rho*rho));
      return -z-std::log(sigma1*sigma2*std::sqrt(1-rho*rho));
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      T mu1=this->parent(0);
      T mu2=this->parent(1);
      T sigma1=this->parent(2);
      T sigma2=this->parent(3);

      if(this->get_current_idx()==0)
	{
	  return std::make_pair(mu1-5*sigma1,mu1+5*sigma1);
	}
      else
	{
	  return std::make_pair(mu2-5*sigma2,mu2+5*sigma2);
	}
    }

    void do_init_value(size_t n) override
    {
      this->set_value(n,this->parent(n));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new bvnormal_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  class bvnormal_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    bvnormal_node_factory()
      :abstract_node_factory<T,T_vector>({"mu1","mu2","sigma1","sigma2","rho"},{"x1","x2"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new bvnormal_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
}

#endif
