#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <memory>
#include <vector>
#include <list>
#include <map>
#include "mcmc_exception.hpp"
#include "stochastic_node.hpp"
#include "deterministic_node.hpp"
#include "observed_node.hpp"
#include "tag.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1,typename T_str>
  class graph
  {
  private:
    std::list<stochastic_node<T_p,T_var1>* > stochastic_node_list;
    std::list<deterministic_node<T_p,T_var1>* > deterministic_node_list;
    std::list<observed_node<T_p,T_var1>* > observed_node_list;
    std::map<tag_t<T_str>,std::shared_ptr<node<T_p,T_var1> > > node_map;

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
      std::vector<T_var1> result;
      int n=0;
      for(auto p=stochastic_node_list.begin();
	  p!=stochastic_node_list.end();++p,++n)
	{
	  for(int i=0;i<(*p)->num_of_dims();++i)
	    {
	      result.push_back((*p)->value(i,0));
	    }
	}
      return result;
    }

    void set_value(const tag_t<T_str>& tag,size_t idx,
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
      ps->set_value(idx,v);
    }

    T_p log_likelihood(const tag_t<T_str>& tag)const
    {
      auto i=node_map.find(tag);
      if(i==node_map.end())
	{
	  throw mcmc_exception("node not found");
	}
      i->second->log_likelihood();
    }

    T_p log_post_prob(const tag_t<T_str>& tag)const
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


    void add_node(node<T_p,T_var1>* pn,
		  const tag_t<T_str>& tag,
		  const std::vector<tag_t<T_str> >& parents,
		  const std::vector<size_t>& idx)
    {
      this->add_node(std::shared_ptr<node<T_p,T_var1> >(pn),tag,parents,idx);
    }

    
    void add_node(const std::shared_ptr<node<T_p,T_var1> >& pn,
		  const tag_t<T_str>& tag,
		  const std::vector<tag_t<T_str> >& parents,
		  const std::vector<size_t>& idx)
    {
      //std::shared_ptr<node<T_p,T_var1> > ptr(pn);
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
      
      auto ps=dynamic_cast<stochastic_node<T_p,T_var1>*>(pn.get());
      auto pd=dynamic_cast<deterministic_node<T_p,T_var1>*>(pn.get());
      auto po=dynamic_cast<observed_node<T_p,T_var1>*>(pn.get());
      if(ps!=nullptr&&pd==nullptr&&po==nullptr)
	{
	  stochastic_node_list.push_back(ps);
	}
      else if(pd!=nullptr&&ps==nullptr&&po==nullptr)
	{
	  deterministic_node_list.push_back(pd);
	}
      else if(po!=nullptr&&pd==nullptr&&ps==nullptr)
	{
	  observed_node_list.push_back(po);
	}
      else
	{
	  throw mcmc_exception("input node is neither stochastic node, nor deterministic node");
	}
      int n=0;
      for(auto& i:parents)
	{
	  auto n_iter=node_map.find(i);
	  pn->connect_to_parent(n_iter->second.get(),n,idx.at(n));
	  ++n;
	}
      node_map[tag]=pn;
    }
  };
}

#endif

