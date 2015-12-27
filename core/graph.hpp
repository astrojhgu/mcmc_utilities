#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <memory>
#include <vector>
#include <list>
#include <map>
#include <stack>
#include <algorithm>
#include <functional>
#include <sstream>
#include "error_handler.hpp"
#include "stochastic_node.hpp"
#include "deterministic_node.hpp"
//#include "tag.hpp"

namespace mcmc_utilities
{
  template <typename T,typename T_tag>
  class graph
  {
  private:
    std::list<stochastic_node<T>* > stochastic_node_list;
    std::list<deterministic_node<T>* > deterministic_node_list;
    std::map<T_tag,std::shared_ptr<node<T> > > node_map;
    std::map<std::shared_ptr<node<T> > ,T_tag, std::owner_less<std::shared_ptr<node<T> > > > tag_map;
  public:
    void copy_from(const graph& rhs)
    {
      std::map<node<T>*,T_tag> tag_map1;
      for(auto i=rhs.tag_map.begin();i!=rhs.tag_map.end();++i)
	{
	  tag_map1[i->first.get()]=i->second;
	}

      std::stack<std::shared_ptr<node<T> > > node_stack;
      std::stack<size_t> parent_id_stack;
      
      for(auto i=rhs.node_map.begin();i!=rhs.node_map.end();++i)
	{
	  if(this->node_map.count(i->first)!=0)
	    {
	      continue;
	    }
	  node_stack.push(i->second);
	  parent_id_stack.push(0);
	  for(;;)
	    {
	      if(parent_id_stack.top()==node_stack.top()->num_of_parents())
		{
		  std::vector<std::pair<T_tag,size_t> > parents(node_stack.top()->num_of_parents());
		  
		  if(node_map.count(rhs.tag_map.find(node_stack.top())->second)==0)
		    {
		      for(size_t i=0;i!=parents.size();++i)
			{
			  parents[i]=std::pair<T_tag,size_t>(tag_map1[node_stack.top()->get_parent(i).first],
							     node_stack.top()->get_parent(i).second);
			}
		      add_node(node_stack.top()->clone(),rhs.tag_map.find(node_stack.top())->second,parents);
		    }
		  node_stack.pop();
		  parent_id_stack.pop();
		  if(node_stack.empty())
		    {
		      break;
		    }
		}
	      else
		{
		  std::shared_ptr<node<T> > p(rhs.node_map.find(tag_map1[node_stack.top()->get_parent(parent_id_stack.top()).first])->second);
		  node_stack.push(p);
		  ++parent_id_stack.top();
		  parent_id_stack.push(0);
		}
	    }
	}
    }
    
    void clear()
    {
      node_map.clear();
      stochastic_node_list.clear();
      deterministic_node_list.clear();
    }

    void sample(const base_urand<T>& rnd)
    {
      for(auto& p:stochastic_node_list)
	{
	  p->sample(rnd);
	}
    }

    void initialize()
    {
      for(auto& p:stochastic_node_list)
	{
	  p->initialize();
	}
    }

    std::vector<T> get_params()const
    {
      std::vector<T> result;
      size_t n=0;
      for(auto p=stochastic_node_list.begin();
	  p!=stochastic_node_list.end();++p,++n)
	{
	  for(size_t i=0;i<(*p)->num_of_dims();++i)
	    {
	      result.push_back((*p)->value(i));
	    }
	}
      return result;
    }

    std::function<T()> get_monitor(const T_tag& tag,size_t n)const
    {
      auto iter=node_map.find(tag);
      if(iter==node_map.end())
	{
	  std::ostringstream oss;
	  oss<<tag;
	  throw node_not_found(oss.str());
	}
      stochastic_node<T>* ps=dynamic_cast<stochastic_node<T>*> (iter->second.get());
      deterministic_node<T>* pd=dynamic_cast<deterministic_node<T>*> (iter->second.get());
      if(ps==nullptr&&pd==nullptr)
	{
	  throw invalid_node_type();
	}
      std::weak_ptr<node<T> > wp(iter->second);
      

      return [n,wp](){
	if(wp.expired())
	  {
	    throw pointer_expired();
	  }
	return wp.lock()->value(n);
      };
    }

    T get_value(const T_tag& tag,size_t idx)
    {
      if(node_map.count(tag)==0)
	{
	  std::ostringstream oss;
	  oss<<tag;
	  throw node_not_found(oss.str());
	}
      stochastic_node<T>* ps=dynamic_cast<stochastic_node<T>* >(node_map[tag].get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      return ps->value(idx);
    }

    void set_value(const T_tag& tag,size_t idx,
		      const T& v)
    {
      if(node_map.count(tag)==0)
	{
	  std::ostringstream oss;
	  oss<<tag;
	  throw node_not_found(oss.str());
	}
      
      stochastic_node<T>* ps=dynamic_cast<stochastic_node<T>* >(node_map[tag].get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      ps->set_value(idx,v);
    }

    T log_likelihood(const T_tag& tag)const
    {
      auto i=node_map.find(tag);
      if(i==node_map.end())
	{
	  std::ostringstream oss;
	  oss<<tag;
	  throw node_not_found(oss.str());
	}
      return i->second->log_likelihood();
    }

    T log_posterior_prob(const T_tag& tag)const
    {
      auto i=node_map.find(tag);
      if(i==node_map.end())
	{
	  std::ostringstream oss;
	  oss<<tag;
	  throw node_not_found(oss.str());
	}

      auto ps=dynamic_cast<stochastic_node<T>*>(i->second.get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      return ps->log_posterior_prob();
    }

    void add_node(node<T>* pn,
		  const T_tag& tag)
    {
      std::vector<std::pair<T_tag,size_t> > parents;
      this->add_node(std::shared_ptr<node<T> >(pn),tag,parents);
    }
    
    void add_node(const std::shared_ptr<node<T> >& pn,
		  const T_tag& tag)
    {
      std::vector<std::pair<T_tag,size_t> > parents;
      this->add_node(pn,tag,parents);
    }

    void add_node(node<T>* pn,
		  const T_tag& tag,
		  const std::vector<std::pair<T_tag,size_t> >& parents
		  )
    {
      this->add_node(std::shared_ptr<node<T> >(pn),tag,parents);
    }
    
    
    void add_node(const std::shared_ptr<node<T> >& pn,
		  const T_tag& tag,
		  const std::vector<std::pair<T_tag,size_t> >& parents
		  )
    {
#if 0
      std::cerr<<"added "<<tag;
      std::cerr<<" with parents:";
      for(auto& i:parents)
	{
	  std::cerr<<"("<<i.first<<":"<<i.second<<") ";
	}
      std::cerr<<std::endl;
#endif 
      //std::shared_ptr<node<T_p,T_var1> > ptr(pn);
      if (node_map.count(tag)!=0)
	{
	  throw node_name_already_used();
	}
      if(tag_map.count(pn)!=0)
	{
	  throw node_already_added();
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
      
      auto ps=dynamic_cast<stochastic_node<T>*>(pn.get());
      auto pd=dynamic_cast<deterministic_node<T>*>(pn.get());
      if(ps!=nullptr&&pd==nullptr)
	{
	  stochastic_node_list.push_back(ps);
	}
      else if(pd!=nullptr&&ps==nullptr)
	{
	  deterministic_node_list.push_back(pd);
	}
      else
	{
	  throw invalid_node_type();
	}
      size_t n=0;
      for(const auto& i:parents)
	{
	  auto n_iter=node_map.find(i.first);
	  pn->connect_to_parent(n_iter->second.get(),n,i.second);
	  ++n;
	}
      node_map[tag]=pn;
      tag_map[pn]=tag;
    }


    void add_node(node<T>*pn,
		  const T_tag& tag,
		  const std::vector<std::pair<std::shared_ptr<node<T> >,size_t> >& parents)
    {
      std::vector<std::pair<node<T>*,size_t> > parents1;
      add_node(std::shared_ptr<node<T> >(pn),tag,parents);
      }
    
    void add_node(const std::shared_ptr<node<T> >& pn,
		  const T_tag& tag,
		  const std::vector<std::pair<std::shared_ptr<node<T> >,size_t> >& parents)
    {
      std::vector<std::pair<T_tag,size_t> > parent_tags;
      for(auto& p:parents)
	{
	  if(tag_map.count(p.first)==0)
	    {
	      throw parents_not_exist();
	    }
	  T_tag tag=tag_map[p.first];
	  parent_tags.push_back({tag,p.second});
	}
      add_node(pn,tag,parent_tags);
    }

#if 0

    void add_node(node<T>* pn,
		  const T_tag& tag,
		  const std::vector<std::pair<node<T>*,size_t> >& parents)
    {
      add_node(std::shared_ptr<node<T> >(pn),tag,parents);
    }

    
    void add_node(const std::shared_ptr<node<T> >& pn,
		  const T_tag& tag,
		  const std::vector<std::pair<node<T>*,size_t> >& parents)
    {
      std::vector<std::pair<T_tag,size_t> > parent_tags;
      for(auto& p:parents)
	{
	  node<T>* pp(p.first);
	  if(tag_map.count(pp)==0)
	    {
	      throw parents_not_exist();
	    }
	  T_tag tag=tag_map[pp];
	  parent_tags.push_back({tag,p.second});
	}
      this->add_node(pn,tag,parent_tags);
    }
#endif
    
  };
}

#endif

