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
#include "step_node.hpp"
#include "abs_node.hpp"
#include "log_node.hpp"
#include "log10_node.hpp"
#include "logit_node.hpp"
#include "ilogit_node.hpp"
#include "sin_node.hpp"
#include "cos_node.hpp"
#include "tan_node.hpp"
#include "arithmetic_node.hpp"
#include "compare_node.hpp"
#include "cond_node.hpp"
#include "phi_node.hpp"
#include "cosmology_nodes.hpp"
#include "ez_node.hpp"
#include <atomic>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class str_node
    :public composed_node<T,std::string,T_vector>
  {
  public:
    static std::map<std::string,std::shared_ptr<abstract_node_factory<T,T_vector> > > node_factories;
    static void init_node_factories()
    {
      std::atomic<bool> initialized(false);

      if(initialized)
	{
	  return;
	}

      //auto& node_factories=str_node::node_factories;
      register_function({"eq"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new eq_node_factory<T,T_vector>()));
      register_function({"lt"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new lt_node_factory<T,T_vector>()));
      register_function({"le"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new le_node_factory<T,T_vector>()));
      register_function({"gt"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new gt_node_factory<T,T_vector>()));
      register_function({"ge"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new ge_node_factory<T,T_vector>()));
      register_function({"add"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new add_node_factory<T,T_vector>()));
      register_function({"sub"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new sub_node_factory<T,T_vector>()));
      register_function({"neg"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new neg_node_factory<T,T_vector>()));
      register_function({"pos"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new pos_node_factory<T,T_vector>()));
      register_function({"mul"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new mul_node_factory<T,T_vector>()));
      register_function({"div"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new div_node_factory<T,T_vector>()));
      register_function({"pow"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new pow_node_factory<T,T_vector>()));
      register_function({"con"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new const_node_factory<T,T_vector>()));
      register_function({"phi"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new phi_node_factory<T,T_vector>()));
      register_function({"sqrt"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new sqrt_node_factory<T,T_vector>()));
      register_function({"step"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new step_node_factory<T,T_vector>()));
      register_function({"abs"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new abs_node_factory<T,T_vector>()));
      register_function({"log"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new log_node_factory<T,T_vector>()));
      register_function({"log10"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new log10_node_factory<T,T_vector>()));
      register_function({"sin"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new sin_node_factory<T,T_vector>()));
      register_function({"cos"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new cos_node_factory<T,T_vector>()));
      register_function({"tan"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new tan_node_factory<T,T_vector>()));
      register_function({"logit"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new logit_node_factory<T,T_vector>()));
      register_function({"ilogit"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new ilogit_node_factory<T,T_vector>()));
      register_function({"D_L"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new luminosity_distance_node_factory<T,T_vector>()));
      register_function({"D_A"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new asize_distance_node_factory<T,T_vector>()));
      register_function({"Ez"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new ez_node_factory<T,T_vector>()));
      register_function({"cond"},std::shared_ptr<abstract_node_factory<T,T_vector> >(new cond_node_factory<T,T_vector>()));
      initialized=true;
    }

    static void register_function(const std::string& name,const std::shared_ptr<abstract_node_factory<T,T_vector> >& ptr)
    {
      node_factories[name]=ptr;
    }

    static void unregister_function(const std::string& name)
    {
      node_factories.erase(name);
    }
  private:
    std::string expression;
    std::vector<std::string> input_names;
  private:
    void count_input(const east::expression_node& en,
		std::set<std::string>& tags)
    {
      for(size_t i=0;i<en.get_num_of_parents();++i)
	{
	  count_input(en.get_parent(i),tags);
	}
      if(en.get_kind()=="var")
	{
	  tags.insert(en.get_symbol());
	}
    }

  private:
    std::shared_ptr<deterministic_node<T,T_vector> > add_node(const east::expression_node& en,int& n)
    {
      if(en.get_kind()=="con")
	{
	  if(this->elements.count(en.get_symbol())==0)
	    {
	      auto p=std::dynamic_pointer_cast<deterministic_node<T,T_vector> >(node_factories["con"]->get_node({(T)std::stod(en.get_symbol())}));
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
	  std::shared_ptr<deterministic_node<T,T_vector> > p;
	  typename std::map<std::string,std::shared_ptr<abstract_node_factory<T,T_vector> > >::iterator iter;
	  if(en.get_kind()!="ftn")
	    {
	      iter=node_factories.find(en.get_kind());
	      //p=std::dynamic_pointer_cast<deterministic_node<T,T_vector> >(node_factories[en.get_kind()]->get_node());
	    }
	  else
	    {
	      iter=node_factories.find(en.get_symbol());
	      //p=std::dynamic_pointer_cast<deterministic_node<T,T_vector> >(node_factories[en.get_symbol()]->get_node());
	    }
	  if(iter==node_factories.end())
	    {
	      throw mcmc_exception("node not registered");
	    }
	  p=std::dynamic_pointer_cast<deterministic_node<T,T_vector> >(iter->second->get_node());
	  if(p==nullptr)
	    {
	      throw mcmc_exception("maybe not a deterministic node");
	    }
	  for(size_t i=0;i<en.get_num_of_parents();++i)
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
    static T eval_expr(const std::string& expression1,
		     const std::vector<std::string>& input_names1,
		     const std::vector<T>& params)
    {
      str_node<T,T_vector> sn(expression1,input_names1);
      std::vector<std::shared_ptr<node<T,T_vector> > > parents;
      for(int i=0;i<input_names1.size();++i)
	{
	  parents.push_back(std::shared_ptr<node<T,T_vector> >(new const_node<T,T_vector>(params[i])));
	  sn.connect_to_parent(parents.back().get(),i,0);
	}
      return sn.value(0); 
    }
		  
  public:
    str_node(const std::string& expression1,
	     const std::vector<std::string>& input_names1)
      :composed_node<T,std::string,T_vector>(input_names1.size(),1),
      expression(expression1),input_names(input_names1)
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
      for(size_t i=0;i<input_names.size();++i)
	{
	  std::shared_ptr<deterministic_node<T,T_vector> > pn(new forward_node<T,T_vector>);
	  this->elements.insert(std::make_pair(input_names[i],pn));
	  this->param_list.push_back(pn);	  	  
	}
      int n=1;
      auto pn=add_node(*enode,n);
      this->return_list.push_back({pn,0});
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      return std::shared_ptr<node<T,T_vector> >(new str_node<T,T_vector>(expression,input_names));
    }
  };

  template <typename T,template <typename TE> class T_vector>
  std::map<std::string,std::shared_ptr<abstract_node_factory<T,T_vector> > > str_node<T,T_vector>::node_factories;
  
  template <typename T,template <typename TE> class T_vector>
  void add_node_to(composed_node<T,std::string,T_vector>& cn,const east::expression_node& en,std::set<std::string>& tag_set)
  {
    
  }
}


#endif
