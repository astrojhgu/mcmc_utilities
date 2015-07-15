#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <memory>
#include <vector>
#include <map>
#include "mcmc_exception.hpp"
#include "gibbs_sampler.hpp"
#include "base_urand.hpp"
#include "stochastic_node.hpp"
#include "deterministic_node.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1,typename T_str>
  class graph
  {
  private:
    std::vector<stochastic_node<T_p,T_var1>* > stochastic_node_list;
    std::vector<deterministic_node<T_p,T_var1>* > deterministic_node_list;
    std::vector<node<T_p,T_var1>* > observed_node_list;
    std::map<std::pair<T_str,int>,std::shared_ptr<node<T_p,T_var1> > > node_map;

  public:
    void clear()
    {
      node_map.clear();
      stochastic_node_list.clear();
      deterministic_node_list.clear();
      observed_node_list.clear();
    }

    void sample(const base_urand<T_p>& urand)
    {
      for(auto& p:stochastic_node_list)
	{
	  p->sample(urand);
	}
    }

    std::vector<T_var1> get_params()const
    {
      std::vector<T_var1> result(stochastic_node_list.size());
      for(int i=0;i<result.size();++i)
	{
	  result[i]=stochastic_node_list[i]->value();
	}
      return result;
    }

    void set_value(const std::pair<T_str,int>& tag,
		      const T_var1& v)
    {
      if(node_map.count(tag)==0)
	{
	  throw mcmc_exception("node not found by tag");
	}
      
      stochastic_node<T_p,T_var1>* ps=dynamic_cast<stochastic_node<T_p,T_var1>* >(node_map[tag].get());
      if(ps==nullptr)
	{
	  throw mcmc_exception("this node is not a stochastic node");
	}
      ps->set_value(v);

    }

    T_p log_likelihood(const std::pair<T_str,int>& tag)const
    {
      auto i=node_map.find(tag);
      if(i==node_map.end())
	{
	  throw mcmc_exception("node not found");
	}
      i->second->log_likelihood();
    }

    T_p log_post_prob(const std::pair<T_str,int>& tag)const
    {
      auto i=node_map.find(tag);
      if(i==node_map.end())
	{
	  throw mcmc_exception("node not found");
	}

      auto ps=dynamic_cast<stochastic_node<T_p,T_var1>*>(i->second.get());
      if(ps==nullptr)
	{
	  throw mcmc_exception("this node is not a stochastic node");
	}
      return ps->log_post_prob();
    }

        
    void add_node(node<T_p,T_var1>* pn,bool observed,
		  const std::pair<T_str,int>& tag,
		  const std::vector<std::pair<T_str,int> >& parents)
    {
      std::shared_ptr<node<T_p,T_var1> > ptr(pn);
      if (node_map.count(tag)!=0)
	{
	  throw mcmc_exception("node name already exists");
	}
      if(pn->num_of_parents()!=parents.size())
	{
	  throw mcmc_exception("parent num do not match");
	}
      for(auto& p:parents)
	{
	  if(node_map.count(p)==0)
	    {
	      throw mcmc_exception("parent should be added before hand");
	    }
	}
      
      auto ps=dynamic_cast<stochastic_node<T_p,T_var1>*>(pn);
      if(ps!=nullptr)
	{
	  if(!observed)
	    {
	      stochastic_node_list.push_back(ps);
	    }
	  else
	    {
	      observed_node_list.push_back(pn);
	    }
	}
      else
	{
	  auto pd=dynamic_cast<deterministic_node<T_p,T_var1>* >(pn);
	  if(pd==nullptr)
	    {
	      throw mcmc_exception("input node is neither stochastic node, nor deterministic node");
	    }
	  deterministic_node_list.push_back(pd);
	}
      int n=0;
      for(auto& i:parents)
	{
	  auto n_iter=node_map.find(i);
	  pn->connect_to_parent(n_iter->second.get(),n);
	  ++n;
	}
      node_map[tag]=ptr;
    }
  };
}

#endif

