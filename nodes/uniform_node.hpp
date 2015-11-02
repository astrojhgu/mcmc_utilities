#ifndef UNIFORM_NODE_HPP
#define UNIFORM_NODE_HPP
#include <core/deterministic_node.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <string>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class uniform_node
    :public stochastic_node<T_p,T_var1>
  {
  private:
    T_var1 a;
    T_var1 b;
  public:
    uniform_node(T_var1 _a,T_var1 _b)
      :stochastic_node<T_p,T_var1>(0,(_a+_b)/2),a(_a),b(_b)
    {}
    
  private:
    T_p do_log_prob()const override
    {
      return 0;
    }
    
    bool is_continuous(size_t)const override
    {
      return true;
    }
    
    std::pair<T_var1,T_var1> do_var_range()const override
    {
      return std::make_pair(a,b);
      //return make_pair(-10,10);
    }
  };

  template <typename T_p,typename T_var1>
  class uniform_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    uniform_node_factory()
      :abstract_node_factory<T_p,T_var1>({},{"x"},{"a","b"},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(
	     const std::vector<T_var1>& scalar_param,
	     const std::vector<std::vector<T_var1> >& vector_param)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new uniform_node<T_p,T_var1>(scalar_param.at(0),scalar_param.at(1)));
    }

    std::string do_get_node_type()const override
    {
      return std::string("stochastic node");
    }
  };
  
  template <typename T_p,typename T_var1>
  class uniform_vnode
    :public vnode<T_p,T_var1>
  {
    T_var1 a;
    T_var1 b;
  public:
    uniform_vnode(std::string n,T_var1 _a,T_var1 _b)
      :vnode<T_p,T_var1>("uniform",n),a(_a),b(_b)
    {
      this->binded=true;
    }
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new uniform_node<T_p,T_var1>(a,b));
    }
    
    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new uniform_vnode<T_p,T_var1>(*this));
    }
  };

  using vuniform=uniform_vnode<double,double>;
};

#endif
