#ifndef UNIFORM_NODE_HPP
#define UNIFORM_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <helper/abstract_node_factory.hpp>
#include <string>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class uniform_node
    :public stochastic_node<T,T_vector>
  {
  private:
    T a;
    T b;
  public:
    uniform_node(T _a,T _b)
      :stochastic_node<T,T_vector>(0,(_a+_b)/2),a(_a),b(_b)
    {}
    
  private:
    T do_log_prob()const override
    {
      return 0;
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T,T> do_var_range()const override
    {
      return std::make_pair(a,b);
      //return make_pair(-10,10);
    }

    void do_init_value(size_t n) override
    {
      this->set_value(0,(a+b)/2);
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      auto p=new uniform_node(a,b);
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }

  };

  template <typename T,template <typename TE> class T_vector>
  class uniform_node_factory
    :public abstract_node_factory<T,T_vector>
  {
  public:
    uniform_node_factory()
      :abstract_node_factory<T,T_vector>({},{"x"},{"a","b"})
    {}
    
  public:
    std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T,T_vector> >(new uniform_node<T,T_vector>(hparam.at(0),hparam.at(1)));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };
  
  template <typename T,template <typename TE> class T_vector>
  class uniform_vnode
    :public vnode<T,T_vector>
  {
    T a;
    T b;
  public:
    uniform_vnode(std::string n,T _a,T _b)
      :vnode<T,T_vector>("uniform",n),a(_a),b(_b)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T,T_vector> > get_node()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new uniform_node<T,T_vector>(a,b));
    }
    
    std::shared_ptr<vnode<T,T_vector> > clone()const override
    {
      return std::shared_ptr<vnode<T,T_vector> >(new uniform_vnode<T,T_vector>(*this));
    }
  };

  //using vuniform=uniform_vnode<double>;
}

#endif
