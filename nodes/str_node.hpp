#ifndef COMPOSE_NODE_HPP
#define COMPOSE_NODE_HPP

#include <memory>
#include <map>
#include <set>
#include <stdexcept>
#include <east.hpp>
#include "core/composed_node.hpp"
#include "phi_node.hpp"
#include "const_node.hpp"
#include "sqrt_node.hpp"
#include "log_node.hpp"
#include "log10_node.hpp"
#include "logit_node.hpp"
#include "ilogit_node.hpp"
#include "sin_node.hpp"
#include "cos_node.hpp"
#include "tan_node.hpp"
#include "arithmetic_node.hpp"
#include "phi_node.hpp"
#include <atomic>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  class str_node
    :public composed_node<T_p,T_var,std::string>
  {
  public:
    static std::map<std::string,std::shared_ptr<abstract_node_factory<T_p,T_var> > > node_factories;
    static void init_node_factories()
    {
      std::atomic<bool> initialized(false);

      if(initialized)
	{
	  return;
	}

      auto& node_factories=str_node::node_factories;

      node_factories["add"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new add_node_factory<T_p,T_var>());
      node_factories["sub"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new sub_node_factory<T_p,T_var>());
      node_factories["neg"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new neg_node_factory<T_p,T_var>());
      node_factories["pos"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new pos_node_factory<T_p,T_var>());
      node_factories["mul"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new mul_node_factory<T_p,T_var>());
      node_factories["div"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new div_node_factory<T_p,T_var>());
      node_factories["pow"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new pow_node_factory<T_p,T_var>());
      node_factories["con"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new const_node_factory<T_p,T_var>());
      node_factories["phi"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new phi_node_factory<T_p,T_var>());
      node_factories["sqrt"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new sqrt_node_factory<T_p,T_var>());
      node_factories["log"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new log_node_factory<T_p,T_var>());
      node_factories["log10"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new log10_node_factory<T_p,T_var>());
      node_factories["sin"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new sin_node_factory<T_p,T_var>());
      node_factories["cos"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new cos_node_factory<T_p,T_var>());
      node_factories["tan"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new tan_node_factory<T_p,T_var>());
      node_factories["logit"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new logit_node_factory<T_p,T_var>());
      node_factories["ilogit"]=std::shared_ptr<abstract_node_factory<T_p,T_var> >(new ilogit_node_factory<T_p,T_var>());
      
      initialized=true;
    }
  private:
    void count_input(const east::expression_node& en,
		std::set<std::string>& tags)
    {
      for(int i=0;i<en.get_num_of_parents();++i)
	{
	  count_input(en.get_parent(i),tags);
	}
      if(en.get_kind()=="var")
	{
	  tags.insert(en.get_symbol());
	}
    }

  private:
    std::shared_ptr<deterministic_node<T_p,T_var> > add_node(const east::expression_node& en,int& n)
    {
      if(en.get_kind()=="con")
	{
	  if(this->elements.count(en.get_symbol())==0)
	    {
	      auto p=std::dynamic_pointer_cast<deterministic_node<T_p,T_var> >(node_factories["con"]->get_node({(T_var)std::stod(en.get_symbol())}));
	      this->elements[en.get_symbol()]=p;
	      return p;
	    }
	  else
	    {
	      return this->elements[en.get_symbol()];
	    }
	}
      else if(en.get_kind()=="var")
	{
	  return this->elements[en.get_symbol()];
	}
      else
	{
	  std::shared_ptr<deterministic_node<T_p,T_var> > p;
	  if(en.get_kind()!="ftn")
	    {
	      p=std::dynamic_pointer_cast<deterministic_node<T_p,T_var> >(node_factories[en.get_kind()]->get_node());
	    }
	  else
	    {
	      p=std::dynamic_pointer_cast<deterministic_node<T_p,T_var> >(node_factories[en.get_symbol()]->get_node());
	    }
	  for(int i=0;i<en.get_num_of_parents();++i)
	    {
	      auto p1=add_node(en.get_parent(i),n);
	      p->connect_to_parent(p1.get(),i,0);
	    }
	  std::string tag="node_"+std::to_string(n++);
	  this->elements[tag]=p;
	  return p;
	}
    }
		  
  public:
    str_node(const std::string& expression,
		    const std::vector<std::string>& input_names)
      :composed_node<T_p,T_var,std::string>(input_names.size(),1)
    {
      str_node::init_node_factories();
      east::parser parser;
      std::shared_ptr<east::expression_node> enode=parser.parse(expression);
      std::set<std::string> tags;
      
      count_input(*enode,tags);
      if(input_names.size()!=tags.size())
	{
	  throw std::logic_error("input names mismatch");
	}
      for(auto s:input_names)
	{
	  if(tags.count(s)!=1)
	    {
	      throw std::logic_error("input names mismatch");
	    }
	}
      for(int i=0;i<input_names.size();++i)
	{
	  std::shared_ptr<deterministic_node<T_p,T_var> > pn(new forward_node<T_p,T_var>);
	  this->elements.insert(std::make_pair(input_names[i],pn));
	  this->param_list.push_back(pn);	  	  
	}
      int n=1;
      auto pn=add_node(*enode,n);
      this->return_list.push_back({pn,0});
    }
  };

  template <typename T_p,typename T_var>
  std::map<std::string,std::shared_ptr<abstract_node_factory<T_p,T_var> > > str_node<T_p,T_var>::node_factories;
  
  template <typename T_p,typename T_var>
  void add_node_to(composed_node<T_p,T_var,std::string>& cn,const east::expression_node& en,std::set<std::string>& tag_set)
  {
    
  }
};


#endif
