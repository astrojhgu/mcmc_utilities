#ifndef BVNORMAL_NODE_HPP
#define BVNORMAL_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T>
  class bvnormal_node
    :public stochastic_node<T>
  {
  private:
  public:
    bvnormal_node()
      :stochastic_node<T>(5,{0,0})
    {}
    
  private:
    T do_log_prob()const override
    {
      const static T PI=std::atan(1.0)*4;
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

    void do_initialize(size_t n) override
    {
      this->set_value(n,this->parent(n));
    }
  };
  
  
  template <typename T>
  class bvnormal_vnode
    :public vnode<T>
  {
  public:
    bvnormal_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T>&,size_t> >& p)
      :vnode<T>("bvnormal",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new bvnormal_node<T>);
    }

    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new bvnormal_vnode<T>(*this));
    }
  };

  template <typename T>
  class bvnormal_node_factory
    :public abstract_node_factory<T>
  {
  public:
    bvnormal_node_factory()
      :abstract_node_factory<T>({"mu1","mu2","sigma1","sigma2","rho"},{"x1","x2"},{})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new bvnormal_node<T>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vbvnormal=bvnormal_vnode<double>;
};

#endif
