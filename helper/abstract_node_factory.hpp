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
    const std::vector<std::string> scalar_hparam_names;
    const std::vector<std::string> vector_hparam_names;
  public:
    abstract_node_factory(const std::vector<std::string>& iname,
			  const std::vector<std::string>& oname,
			  const std::vector<std::string>& sname,
			  const std::vector<std::string>& vname)
      :input_names(iname),
       output_names(oname),
       scalar_hparam_names(sname),
       vector_hparam_names(vname)
    {}

    std::shared_ptr<node<T_p,T_var1> >
    get_node(
	     const std::vector<T_var1>& scalar_param,
	     const std::vector<std::vector<T_var1> >& vector_param)const
    {
      if(scalar_param.size()!=scalar_hparam_names.size())
	{
	  throw mcmc_exception("scalar param number mismatch");
	}
      
      if(vector_param.size()!=vector_hparam_names.size())
	{
	  throw mcmc_exception("vector param number mismatch");
	}
      return do_get_node(scalar_param,vector_param);
    }

    std::shared_ptr<node<T_p,T_var1> >
    get_node(const std::vector<T_var1>& scalar_param)const
    {
      return get_node(scalar_param,{});
    }

    std::shared_ptr<node<T_p,T_var1> >
    get_node(const std::vector<std::vector<T_var1> >& vector_param)const
    {
      return get_node({},vector_param);
    }

    std::shared_ptr<node<T_p,T_var1> >
    get_node()const
    {
      return get_node({},{});
    }
    
  private:
    virtual std::shared_ptr<node<T_p,T_var1> >
    do_get_node(
		const std::vector<T_var1>& scalar_param,
		const std::vector<std::vector<T_var1> >& vector_param)const=0;
  };
}


#endif
