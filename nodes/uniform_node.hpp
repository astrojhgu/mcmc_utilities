#ifndef UNIFORM_NODE_HPP
#define UNIFORM_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <helper/abstract_node_factory.hpp>
#include <string>

namespace mcmc_utilities
{
  template <typename T>
  class uniform_node
    :public stochastic_node<T>
  {
  private:
    T a;
    T b;
  public:
    uniform_node(T _a,T _b)
      :stochastic_node<T>(0,(_a+_b)/2),a(_a),b(_b)
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

    void do_initialize(size_t n)override
    {
      this->set_value(0,(a+b)/2);
    }
  };

  template <typename T>
  class uniform_node_factory
    :public abstract_node_factory<T>
  {
  public:
    uniform_node_factory()
      :abstract_node_factory<T>({},{"x"},{"a","b"})
    {}
    
  public:
    std::shared_ptr<node<T> >
    do_get_node(const std::vector<T>& hparam)const override
    {
      return std::shared_ptr<node<T> >(new uniform_node<T>(hparam.at(0),hparam.at(1)));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };
  
  template <typename T>
  class uniform_vnode
    :public vnode<T>
  {
    T a;
    T b;
  public:
    uniform_vnode(std::string n,T _a,T _b)
      :vnode<T>("uniform",n),a(_a),b(_b)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T> > get_node()const override
    {
      return std::shared_ptr<node<T> >(new uniform_node<T>(a,b));
    }
    
    std::shared_ptr<vnode<T> > clone()const override
    {
      return std::shared_ptr<vnode<T> >(new uniform_vnode<T>(*this));
    }
  };

  using vuniform=uniform_vnode<double>;
};

#endif
