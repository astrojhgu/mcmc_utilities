#ifndef BVNORMAL_NODE_HPP
#define BVNORMAL_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class bvnormal_node
    :public stochastic_node<T_p,T_var1>
  {
  private:
  public:
    bvnormal_node()
      :stochastic_node<T_p,T_var1>(5,{0,0})
    {}
    
  private:
    T_p do_log_prob()const override
    {
      const static T_var1 PI=std::atan(1.0)*4;
      T_var1 x1=this->value(0,0);
      T_var1 x2=this->value(1,0);
      T_var1 mu1=this->parent(0,0);
      T_var1 mu2=this->parent(1,0);
      T_var1 sigma1=this->parent(2,0);
      T_var1 sigma2=this->parent(3,0);
      T_var1 rho=this->parent(4,0);
      //return -(x-mu)*(x-mu)/(2*sigma*sigma)-std::log(sigma*std::sqrt(2*PI));
      T_var1 X1=(x1-mu1)/sigma1;
      T_var1 X2=(x2-mu2)/sigma2;
      T_var1 z=(X1*X1+X2*X2-2*rho*X1*X2)/(2*(1-rho*rho));
      return -z-std::log(sigma1*sigma2*std::sqrt(1-rho*rho));
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T_var1,T_var1> do_var_range()const override
    {
      T_var1 mu1=this->parent(0,0);
      T_var1 mu2=this->parent(1,0);
      T_var1 sigma1=this->parent(2,0);
      T_var1 sigma2=this->parent(3,0);

      if(this->get_current_idx()==0)
	{
	  return std::make_pair(mu1-5*sigma1,mu1+5*sigma1);
	}
      else
	{
	  return std::make_pair(mu2-5*sigma2,mu2+5*sigma2);
	}
    }
  };
  
  
  template <typename T_p,typename T_var1>
  class bvnormal_vnode
    :public vnode<T_p,T_var1>
  {
  public:
    bvnormal_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T_p,T_var1>&,size_t> >& p)
      :vnode<T_p,T_var1>("bvnormal",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new bvnormal_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new bvnormal_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  class bvnormal_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    bvnormal_node_factory()
      :abstract_node_factory<T_p,T_var1>({"mu1","mu2","sigma1","sigma2","rho"},{"x1","x2"},{},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(
	     const std::vector<T_var1>& scalar_param,
	     const std::vector<std::vector<T_var1> >& vector_param)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new bvnormal_node<T_p,T_var1>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vbvnormal=bvnormal_vnode<double,double>;
};

#endif
