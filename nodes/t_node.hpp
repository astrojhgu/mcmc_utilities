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
  template <typename T,template <typename TE> class T_vector>
  class t_node
    :public stochastic_node<T,T_vector>
  {
  private:
  public:
    t_node()
      :stochastic_node<T,T_vector>(3,0)
    {}
    
  private:
    T do_log_prob()const override
    {
      const static T PI=std::atan(1.0)*4;
      T x=this->value(0);
      T mu=this->parent(0);
      T sigma=this->parent(1);
      T tau=1/(sigma*sigma);
      T k=this->parent(2);
      //return -(x-mu)*(x-mu)/(2*sigma*sigma)-std::log(sigma*std::sqrt(2*PI));
      return std::lgamma((k+1)/2)-std::lgamma(k/2)+std::log(tau/k/PI)/2-
	(k+1)/2*std::log(1+tau*(x-mu)*(x-mu)/k);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      T mu=this->parent(0);
      T sigma=std::abs(this->parent(1));

      //return std::make_pair(mu-5*sigma,mu+5*sigma);
      //return make_pair(-10,10);
      T n=0;
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

    void do_initialize(size_t n)override
    {
      this->set_value(0,this->parent(0));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new t_node;
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };
  
  
  template <typename T,template <typename TE> class T_vector>
  class t_vnode
    :public vnode<T,T_vector>
  {
  public:
    t_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T,T_vector>&,size_t> >& p)
      :vnode<T,T_vector>("t",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new t_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new t_vnode<T,T_vector>(*this));
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class t_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    t_node_factory()
      :abstract_node_factory<T,T_vector>({"mu","sigma","k"},{"x"},{})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new t_node<T,T_vector>);
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  //using vt=t_vnode<double>;
};

#endif
