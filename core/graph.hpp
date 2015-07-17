#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <memory>
#include <vector>
#include <map>
#include "mcmc_exception.hpp"
#include "stochastic_node.hpp"
#include "deterministic_node.hpp"

namespace mcmc_utilities
{
  template <typename T_str>
  class tag_t
  {
  private:
    T_str primary_name;
    std::vector<int> index;
  public:
    tag_t()=delete;
    tag_t(const T_str& n,const std::vector<int>& idx)
      :primary_name(n),index(idx)
    {}

    tag_t(const T_str& n)
      :primary_name(n),index()
    {
      
    }

  public:
    bool operator<(const tag_t<T_str>& rhs)const
    {
      //return primary_name>=rhs.primary_name?true:
      if(primary_name==rhs.primary_name)
	{
	  if(index.size()!=rhs.index.size())
	    {
	      return index.size()<rhs.index.size();
	    }
	  for(int i=0;i<index.size();++i)
	    {
	      if(index[i]==rhs.index[i])
		{
		  continue;
		}
	      return index[i]<rhs.index[i];
	    }
	}
      return (primary_name<rhs.primary_name);
    }
  };

  template <T_str>
  tag_t form_tag(const T_str& str)
  {
    
  }
  
  
  template <typename T_p,typename T_var1,typename T_str>
  class graph
  {
  private:
    std::vector<stochastic_node<T_p,T_var1>* > stochastic_node_list;
    std::vector<deterministic_node<T_p,T_var1>* > deterministic_node_list;
    std::vector<node<T_p,T_var1>* > observed_node_list;
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
      std::vector<T_var1> result(stochastic_node_list.size());
      for(int i=0;i<result.size();++i)
	{
	  result[i]=stochastic_node_list[i]->value();
	}
      return result;
    }

    void set_value(const tag_t<T_str>& tag,
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


    void add_node(node<T_p,T_var1>* pn,bool observed,
		  const tag_t<T_str>& tag,
		  const std::vector<tag_t<T_str> >& parents)
    {
      this->add_node(std::shared_ptr<node<T_p,T_var1> >(pn),observed,tag,parents);
    }

    
    void add_node(const std::shared_ptr<node<T_p,T_var1> >& pn,bool observed,
		  const tag_t<T_str>& tag,
		  const std::vector<tag_t<T_str> >& parents)
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
      if(ps!=nullptr)
	{
	  if(!observed)
	    {
	      stochastic_node_list.push_back(ps);
	    }
	  else
	    {
	      observed_node_list.push_back(pn.get());
	    }
	}
      else
	{
	  auto pd=dynamic_cast<deterministic_node<T_p,T_var1>* >(pn.get());
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
      node_map[tag]=pn;
    }
  };
}

#endif

