#ifndef TRUNC_PARETO_NODE_HPP
#define TRUNC_PARETO_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>
#include <helper/abstract_node_factory.hpp>
#include <math/distributions.hpp>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class trunc_pareto_node
    :public stochastic_node<T,T_vector>
  {
  private:
    T xmax;
  public:
    trunc_pareto_node(T xm)
      :stochastic_node<T,T_vector>(2,0),xmax(xm)
    {}
    
  private:
    T do_log_prob()const override
    {
      T x=this->value(0);
      T c=this->parent(0);
      T alpha=this->parent(1);
      return logdpar<T,T>(x,c,alpha);
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      T c=this->parent(0);
      T alpha=this->parent(1);

      T eta=1e-4;
     
      return std::make_pair(c,std::min(c*std::pow(eta,-(T)(1.)/alpha),xmax));
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,this->parent(0));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new trunc_pareto_node(xmax);
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };
  
  
  template <typename T,template <typename TE> class T_vector>
  class trunc_pareto_vnode
    :public vnode<T,T_vector>
  {
  public:
    trunc_pareto_vnode(std::string n,const std::initializer_list<std::pair<const vnode<T,T_vector>&,size_t> >& p)
      :vnode<T,T_vector>("trunc_pareto",n,p)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new trunc_pareto_node<T,T_vector>);
    }

    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new trunc_pareto_vnode<T,T_vector>(*this));
    }
  };

  template <typename T,template <typename TE> class T_vector>
  class trunc_pareto_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    trunc_pareto_node_factory()
      :abstract_node_factory<T,T_vector>({"mu","sigma"},{"x"},{"xmax"})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new trunc_pareto_node<T,T_vector>(hparam[0]));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }

  };
  
  //using vtrunc_pareto=trunc_pareto_vnode<double>;
};

#endif
