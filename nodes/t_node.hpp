#ifndef T_NODE_HPP
#define T_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <cmath>
#include <helper/abstract_node_factory.hpp>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class t_node
    :public stochastic_node<T_p,T_var1>
  {
  private:
  public:
    t_node()
      :stochastic_node<T_p,T_var1>(3,0)
    {}
    
  private:
    T_p do_log_prob()const override
    {
      const static T_var1 PI=std::atan(1.0)*4;
      T_var1 x=this->value(0);
      T_var1 mu=this->parent(0);
      T_var1 sigma=this->parent(1);
      T_var1 tau=1/(sigma*sigma);
      T_var1 k=this->parent(2);
      //return -(x-mu)*(x-mu)/(2*sigma*sigma)-std::log(sigma*std::sqrt(2*PI));
      return std::lgamma((k+1)/2)-std::lgamma(k/2)+std::log(tau/k/PI)/2-
	(k+1)/2*std::log(1+tau*(x-mu)*(x-mu)/k);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T_var1,T_var1> do_var_range()const override
    {
      T_var1 mu=this->parent(0);
      T_var1 sigma=std::abs(this->parent(1));

      //return std::make_pair(mu-5*sigma,mu+5*sigma);
      //return make_pair(-10,10);
      T_var1 n=0;
      switch((int)this->parent(2))
	{
	case 1:
	  n=500;
	  break;
	case 2:
	  n=30;
	  break;
	case 3:
	  n=10;
	  break;
	default:
	  n=5;
	  break;
	};
      return std::make_pair(mu-n*sigma,mu+n*sigma);
    }
  };
  
  
  template <typename T_p,typename T_var1>
  class t_vnode
    :public vnode<T_p,T_var1>
  {
  public:
    t_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T_p,T_var1>&,size_t> >& p)
      :vnode<T_p,T_var1>("t",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new t_node<T_p,T_var1>);
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new t_vnode<T_p,T_var1>(*this));
    }
  };

  template <typename T_p,typename T_var1>
  class t_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    t_node_factory()
      :abstract_node_factory<T_p,T_var1>({"mu","sigma","k"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(const std::vector<T_var1>& hparam)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new t_node<T_p,T_var1>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  using vt=t_vnode<double,double>;
};

#endif
