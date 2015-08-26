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
//#include "tag.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1,typename T_tag>
  class graph
  {
  private:
    std::list<stochastic_node<T_p,T_var1>* > stochastic_node_list;
    std::list<deterministic_node<T_p,T_var1>* > deterministic_node_list;
    std::list<observed_node<T_p,T_var1>* > observed_node_list;
    std::map<T_tag,std::shared_ptr<node<T_p,T_var1> > > node_map;

  public:
    void clear()
    {
      node_map.clear();
      stochastic_node_list.clear();
      deterministic_node_list.clear();
      observed_node_list.clear();
    }

    void sample(const base_urand<T_p>& rnd)
    {
      for(auto& p:stochastic_node_list)
	{
	  p->sample(rnd);
	}
    }

    std::vector<T_var1> get_params()const
    {
      std::vector<T_var1> result;
      size_t n=0;
      for(auto p=stochastic_node_list.begin();
	  p!=stochastic_node_list.end();++p,++n)
	{
	  for(size_t i=0;i<(*p)->num_of_dims();++i)
	    {
	      result.push_back((*p)->value(i,0));
	    }
	}
      return result;
    }

    auto get_monitor(const T_tag& tag,size_t n)const
    {
      auto iter=node_map.find(tag);
      if(iter==node_map.end())
	{
	  throw node_not_found();
	}
      stochastic_node<T_p,T_var1>* ps=dynamic_cast<stochastic_node<T_p,T_var1>*> (iter->second.get());
      deterministic_node<T_p,T_var1>* pd=dynamic_cast<deterministic_node<T_p,T_var1>*> (iter->second.get());
      if(ps==nullptr&&pd==nullptr)
	{
	  throw invalid_node_type();
	}
      std::weak_ptr<node<T_p,T_var1> > wp(iter->second);
      

      return [n,wp](){
	if(wp.expired())
	  {
	    throw pointer_expired();
	  }
	return wp.lock()->value(n,0);
      };
    }

    T_var1 get_value(const T_tag& tag,size_t idx)
    {
      if(node_map.count(tag)==0)
	{
	  throw node_not_found();
	}
      stochastic_node<T_p,T_var1>* ps=dynamic_cast<stochastic_node<T_p,T_var1>* >(node_map[tag].get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      return ps->value(idx,0);
    }

    void set_value(const T_tag& tag,size_t idx,
		      const T_var1& v)
    {
      if(node_map.count(tag)==0)
	{
	  throw node_not_found();
	}
      
      stochastic_node<T_p,T_var1>* ps=dynamic_cast<stochastic_node<T_p,T_var1>* >(node_map[tag].get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      ps->set_value(idx,v);
    }

    T_p log_likelihood(const T_tag& tag)const
    {
      auto i=node_map.find(tag);
      if(i==node_map.end())
	{
	  throw node_not_found();
	}
      i->second->log_likelihood();
    }

    T_p log_post_prob(const T_tag& tag)const
    {
      auto i=node_map.find(tag);
      if(i==node_map.end())
	{
	  throw node_not_found();
	}

      auto ps=dynamic_cast<stochastic_node<T_p,T_var1>*>(i->second.get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      return ps->log_post_prob();
    }


    void add_node(node<T_p,T_var1>* pn,
		  const T_tag& tag,
		  const std::vector<std::pair<T_tag,size_t> >& parents
		  )
    {
      this->add_node(std::shared_ptr<node<T_p,T_var1> >(pn),tag,parents);
    }
    
    
    void add_node(const std::shared_ptr<node<T_p,T_var1> >& pn,
		  const T_tag& tag,
		  const std::vector<std::pair<T_tag,size_t> >& parents
		  )
    {
      //std::shared_ptr<node<T_p,T_var1> > ptr(pn);
      if (node_map.count(tag)!=0)
	{
	  throw node_name_already_used();
	}
      if(pn->num_of_parents()!=parents.size())
	{
	  throw parent_num_mismatch();
	}
      for(auto& p:parents)
	{
	  if(node_map.count(p.first)==0)
	    {
	      throw parents_not_exist();
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
	  throw invalid_node_type();
	}
      size_t n=0;
      for(auto& i:parents)
	{
	  auto n_iter=node_map.find(i.first);
	  pn->connect_to_parent(n_iter->second.get(),n,i.second);
	  ++n;
	}
      node_map[tag]=pn;
    }
  };
}

#endif

