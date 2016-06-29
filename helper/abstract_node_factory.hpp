#ifndef ABSTRACT_NODE_FACTORY
#define ABSTRACT_NODE_FACTORY

#include <core/node.hpp>
#include <memory>
#include <string>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class abstract_node_factory
  {
  public:
    const T_vector<std::string> input_names;
    const T_vector<std::string> output_names;
    const T_vector<std::string> hparam_names;
    
  public:
    abstract_node_factory(const T_vector<std::string>& iname,
			  const T_vector<std::string>& oname,
			  const T_vector<std::string>& hname)
      :input_names(iname),
       output_names(oname),
       hparam_names(hname)
    {}

    virtual ~abstract_node_factory(){}

    std::shared_ptr<node<T,T_vector> >
    get_node(const T_vector<T>& hparam)const
    {
      if(get_size(hparam)!=get_size(hparam_names))
	{
	  throw mcmc_exception("param number mismatch");
	}
      
      return do_get_node(hparam);
    }

    
    std::shared_ptr<node<T,T_vector> >
    get_node()const
    {
      return get_node({});
    }

  public:
    std::string get_node_type()const
    {
      return do_get_node_type();
    }

    T_vector<std::string> get_hparam_names()const
    {
      return hparam_names;
    }
    
  private:
    virtual std::shared_ptr<node<T,T_vector> >
    do_get_node(const T_vector<T>& hparam)const=0;

    virtual std::string do_get_node_type()const=0;
    
  };

  
  template <typename T,template <typename TE> class T_vector>
  std::shared_ptr<stochastic_node<T,T_vector> > to_stochastic(const std::shared_ptr<node<T,T_vector> >& sp)
  {
    return std::dynamic_pointer_cast<stochastic_node<T,T_vector> >(sp);
  }

  template <typename T,template <typename TE> class T_vector>
  std::shared_ptr<deterministic_node<T,T_vector> > to_deterministic(const std::shared_ptr<node<T,T_vector> >& sp)
  {
    return std::dynamic_pointer_cast<deterministic_node<T,T_vector> >(sp);
  }
}


#endif
