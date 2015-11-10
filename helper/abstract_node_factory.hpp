#ifndef ABSTRACT_NODE_FACTORY
#define ABSTRACT_NODE_FACTORY

#include <core/node.hpp>
#include <memory>
#include <string>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class abstract_node_factory
  {
  public:
    const std::vector<std::string> input_names;
    const std::vector<std::string> output_names;
    const std::vector<std::string> hparam_names;
    
  public:
    abstract_node_factory(const std::vector<std::string>& iname,
			  const std::vector<std::string>& oname,
			  const std::vector<std::string>& hname)
      :input_names(iname),
       output_names(oname),
       hparam_names(hname)
    {}

    std::shared_ptr<node<T_p,T_var1> >
    get_node(const std::vector<T_var1>& hparam)const
    {
      if(hparam.size()!=hparam_names.size())
	{
	  throw mcmc_exception("param number mismatch");
	}
      
      return do_get_node(hparam);
    }

    
    std::shared_ptr<node<T_p,T_var1> >
    get_node()const
    {
      return get_node({});
    }

  public:
    std::string get_node_type()const
    {
      return do_get_node_type();
    }

    std::vector<std::string> get_hparam_names()const
    {
      return hparam_names;
    }
    
  private:
    virtual std::shared_ptr<node<T_p,T_var1> >
    do_get_node(const std::vector<T_var1>& hparam)const=0;

    virtual std::string do_get_node_type()const=0;
    
  };
}


#endif
