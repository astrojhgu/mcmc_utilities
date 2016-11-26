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
  template <typename T,typename T_tag,template <typename TE> class T_vector>
  class graph
  {
  private:
    std::list<stochastic_node<T,T_vector>* > stochastic_node_list;
    std::list<deterministic_node<T,T_vector>* > deterministic_node_list;
    std::map<T_tag,std::shared_ptr<node<T,T_vector> > > node_map;
    //std::map<std::shared_ptr<node<T,T_vector> > ,T_tag, std::owner_less<std::shared_ptr<node<T,T_vector> > > > tag_map;
    std::map<std::shared_ptr<node<T,T_vector> > ,T_tag > tag_map;
    bool topology_frozen;
    int verbose_level;
  public:
    graph()
      :stochastic_node_list(),
       deterministic_node_list(),
       node_map(),
       tag_map(),
       topology_frozen(false),verbose_level(0)
    {}

    graph(graph&& rhs)
      :stochastic_node_list(rhs.stochastic_node_list),
       deterministic_node_list(rhs.stochastic_node_list),
       node_map(rhs.node_map),
       tag_map(rhs.tag_map),
       topology_frozen(rhs.topology_frozen),
       verbose_level(rhs.verbose_level)
    {
    }

    virtual ~graph(){}

    graph(const graph<T,T_tag,T_vector>&)=delete;
    graph<T,T_tag,T_vector>& operator=(const graph<T,T_tag,T_vector>&)=delete;
    
  public:
    void copy_from(const graph& rhs)
    {
      std::map<node<T,T_vector>*,T_tag> tag_map1;
      for(auto i=std::begin(rhs.tag_map);i!=std::end(rhs.tag_map);++i)
	{
	  tag_map1[i->first.get()]=i->second;
	}

      std::stack<std::shared_ptr<node<T,T_vector> > > node_stack;
      std::stack<size_t> parent_id_stack;
      
      for(auto i=std::begin(rhs.node_map);i!=std::end(rhs.node_map);++i)
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
		  T_vector<std::pair<T_tag,size_t> > parents(node_stack.top()->num_of_parents());
		  
		  if(node_map.count(rhs.tag_map.find(node_stack.top())->second)==0)
		    {
		      for(size_t i=0;i!=get_size(parents);++i)
			{
			  set_element(parents,i,std::pair<T_tag,size_t>(tag_map1[node_stack.top()->get_parent(i).first],
									node_stack.top()->get_parent(i).second));
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
		  std::shared_ptr<node<T,T_vector> > p(rhs.node_map.find(tag_map1[node_stack.top()->get_parent(parent_id_stack.top()).first])->second);
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

    void set_verbose_level(int n)
    {
      verbose_level=n;
    }

    void freeze_topology()
    {
      for(auto& p:stochastic_node_list)
	{
	  p->freeze_topology();
	}
      for(auto& p:deterministic_node_list)
	{
	  p->freeze_topology();
	}
      topology_frozen=true;
    }

    void topology_changed()
    {
      topology_frozen=false;
    }

    void sample(base_urand<T>& rnd)
    {
      if(!topology_frozen)
	{
	  freeze_topology();
	}
      
      stochastic_node<T,T_vector>* p_current=nullptr;
      int n=0;
      //T_vector<stochastic_node<T,T_vector>*> stochastic_node_vector;
      //reserve(stochastic_node_vector,get_size(stochastic_node_list));
      //std::for_each(std::begin(stochastic_node_list),
      //std::end(stochastic_node_list),
      //[&](stochastic_node<T,T_vector>* p){push_back(stochastic_node_vector,p);}
      //);
      try
	{
	  for(auto& p:stochastic_node_list)
	    {
	      if(p->num_of_unobserved()>0)
		{
		  if(verbose_level>=1)
		    {		  
		      std::cerr<<"sampling "<<n+1<<"-th node "<<get_tag(p)<<std::endl;
		    }
		  p_current=p;
		  p->sample(rnd);
		  ++n;
		}
	    }
	}
      catch(mcmc_exception& e)
	{
	  std::ostringstream oss;
	  oss<<"##################\n";
	  oss<<"When sampling\n";
	  oss<<this->get_tag(p_current);
	  oss<<"\n";
	  oss<<__FILE__<<":"<<__LINE__<<"\n";
	  oss<<"##################\n";
	  e.attach_message(oss.str());
	  throw e;
	}
    }

    T log_joint_prob()const
    {
      T result=static_cast<T>(0);
      for(auto& p:stochastic_node_list)
	{
	  result+=p->log_prob();
	}
      return result;
    }

    void initialize()
    {
      for(auto& p:stochastic_node_list)
	{
	  p->init_value();
	}
    }

    T_vector<T> get_params()const
    {
      T_vector<T> result;
      size_t n=0;
      for(auto p=std::begin(stochastic_node_list);
	  p!=std::end(stochastic_node_list);++p,++n)
	{
	  for(size_t i=0;i<(*p)->num_of_dims();++i)
	    {
	      if(!(*p)->is_observed(i))
		{
		  push_back(result,(*p)->value(i));
		}
	    }
	}
      return result;
    }

    void set_params(const T_vector<T>& param)
    {
      size_t n=0;
      for(auto p=std::begin(stochastic_node_list);
	  p!=std::end(stochastic_node_list);++p)
	{
	  for(size_t i=0;i<(*p)->num_of_dims();++i)
	    {
	      if(!(*p)->is_observed(i))
		{
		  (*p)->set_value(i,get_element(param,n++));
		}
	    }
	}
    }

    T eval_logprob(const T_vector<T>& param)
    {
      T result=0;
      size_t n=0;
      for(auto p=std::begin(stochastic_node_list);
	  p!=std::end(stochastic_node_list);++p)
	{
	  for(size_t i=0;i<(*p)->num_of_dims();++i)
	    {
	      if(!(*p)->is_observed(i))
		{
		  T x=get_element(param,n++);
		  (*p)->set_current_idx(i);
		  auto xrange=(*p)->var_range();
		  if(x<=xrange.first||x>=xrange.second)
		    {
		      return -std::numeric_limits<T>::infinity();
		    }
		  (*p)->set_value(i,x);
		}
	    }
	}
      return this->log_joint_prob();
    }
    
    std::function<T()> get_monitor(const T_tag& tag,size_t n)const
    {
      auto iter=node_map.find(tag);
      if(iter==std::end(node_map))
	{
	  std::ostringstream oss;
	  oss<<tag;
	  throw node_not_found(oss.str());
	}
      stochastic_node<T,T_vector>* ps=dynamic_cast<stochastic_node<T,T_vector>*> (iter->second.get());
      deterministic_node<T,T_vector>* pd=dynamic_cast<deterministic_node<T,T_vector>*> (iter->second.get());
      if(ps==nullptr&&pd==nullptr)
	{
	  throw invalid_node_type();
	}
      std::weak_ptr<node<T,T_vector> > wp(iter->second);
      

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
      stochastic_node<T,T_vector>* ps=dynamic_cast<stochastic_node<T,T_vector>* >(node_map[tag].get());
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
      
      stochastic_node<T,T_vector>* ps=dynamic_cast<stochastic_node<T,T_vector>* >(node_map[tag].get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      ps->set_value(idx,v);
    }

    bool is_observed(const T_tag& tag,size_t idx)
    {
      if(node_map.count(tag)==0)
	{
	  std::ostringstream oss;
	  oss<<tag;
	  throw node_not_found(oss.str());
	}
      stochastic_node<T,T_vector>* ps=dynamic_cast<stochastic_node<T,T_vector>* >(node_map[tag].get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      return ps->is_observed(idx);
    }

    void set_observed(const T_tag& tag,size_t idx,bool b)
    {
      if(node_map.count(tag)==0)
	{
	  std::ostringstream oss;
	  oss<<tag;
	  throw node_not_found(oss.str());
	}
      
      stochastic_node<T,T_vector>* ps=dynamic_cast<stochastic_node<T,T_vector>* >(node_map[tag].get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      ps->set_observed(idx,b);
    }

    T log_likelihood(const T_tag& tag)const
    {
      auto i=node_map.find(tag);
      if(i==std::end(node_map))
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
      if(i==std::end(node_map))
	{
	  std::ostringstream oss;
	  oss<<tag;
	  throw node_not_found(oss.str());
	}

      auto ps=dynamic_cast<stochastic_node<T,T_vector>*>(i->second.get());
      if(ps==nullptr)
	{
	  throw invalid_node_type();
	}
      return ps->log_posterior_prob();
    }

    void add_node(node<T,T_vector>* pn,
		  const T_tag& tag)
    {
      T_vector<std::pair<T_tag,size_t> > parents;
      this->add_node(std::shared_ptr<node<T,T_vector> >(pn),tag,parents);
    }
    
    void add_node(const std::shared_ptr<node<T,T_vector> >& pn,
		  const T_tag& tag)
    {
      T_vector<std::pair<T_tag,size_t> > parents;
      this->add_node(pn,tag,parents);
    }

    void add_node(node<T,T_vector>* pn,
		  const T_tag& tag,
		  const T_vector<std::pair<T_tag,size_t> >& parents
		  )
    {
      this->add_node(std::shared_ptr<node<T,T_vector> >(pn),tag,parents);
    }
    
    
    void add_node(const std::shared_ptr<node<T,T_vector> >& pn,
		  const T_tag& tag,
		  const T_vector<std::pair<T_tag,size_t> >& parents
		  )
    {
      topology_changed();
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
      if(pn->num_of_parents()!=get_size(parents))
	{
	  throw parent_num_mismatch();
	}
      for(auto& p:parents)
	{
	  if(node_map.count(p.first)==0)
	    {
	      std::cerr<<"parent to be added:"<<p.first<<std::endl;
	      throw parents_not_exist();
	    }
	}
      
      auto ps=dynamic_cast<stochastic_node<T,T_vector>*>(pn.get());
      auto pd=dynamic_cast<deterministic_node<T,T_vector>*>(pn.get());
      if(ps!=nullptr&&pd==nullptr)
	{
	  push_back(stochastic_node_list,ps);
	}
      else if(pd!=nullptr&&ps==nullptr)
	{
	  push_back(deterministic_node_list,pd);
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


    void add_node(node<T,T_vector>*pn,
		  const T_tag& tag,
		  const T_vector<std::pair<std::shared_ptr<node<T,T_vector> >,size_t> >& parents)
    {
      add_node(std::shared_ptr<node<T,T_vector> >(pn),tag,parents);
    }
    
    void add_node(const std::shared_ptr<node<T,T_vector> >& pn,
		  const T_tag& tag,
		  const T_vector<std::pair<std::shared_ptr<node<T,T_vector> >,size_t> >& parents)
    {
      T_vector<std::pair<T_tag,size_t> > parent_tags;
      for(auto& p:parents)
	{
	  if(tag_map.count(p.first)==0)
	    {
	      throw parents_not_exist();
	    }
	  T_tag tag=tag_map[p.first];
	  push_back(parent_tags,std::pair<T_tag,size_t>{tag,p.second});
	}
      add_node(pn,tag,parent_tags);
    }

    T_tag get_tag(const std::shared_ptr<node<T,T_vector> >& p)const
    {
      auto result= tag_map.find(p);
      if(result==std::end(tag_map))
	{
	  throw node_not_found();
	}
      return result->second;
    }

    T_tag get_tag(const node<T,T_vector>* p)const
    {
      return get_tag(std::shared_ptr<node<T,T_vector> >(const_cast<node<T,T_vector>*>(p),[](node<T,T_vector>*){}));
    }

    std::map<T_tag,T_vector<std::pair<T_tag,size_t> > > topology()const
    {
      std::map<T_tag,T_vector<std::pair<T_tag,size_t> > >  result;
      for(auto& p:node_map)
	{
	  auto tag=p.first;
	  auto pnode=p.second;
	  size_t nparents=pnode->num_of_parents();
	  T_vector<std::pair<T_tag,size_t> > parent_list;
	  for(size_t i=0;i<nparents;++i)
	    {
	      auto px=pnode->get_parent(i);
	      T_tag ptag=get_tag(px.first);
	      push_back(parent_list,std::pair<T_tag,size_t>{ptag,px.second});
	    }
	  result[tag]=parent_list;
	}
      return result;
    }

    //const node<T,T_vector>* get_node(const T_tag& t)const
    std::weak_ptr<node<T,T_vector> > get_node(const T_tag& t)const
    {
      auto i=node_map.find(t);
      if(i==std::end(node_map))
	{
	  std::ostringstream oss;
	  oss<<t;
	  throw node_not_found(oss.str());
	}
      return std::weak_ptr<node<T,T_vector> >(i->second);
    }
  };

  template <typename T,typename T_tag,template <typename TE> class T_vector>
  graph<T,T_tag,T_vector> clone(const graph<T,T_tag,T_vector>& rhs )
  {
    graph<T,T_tag,T_vector> g;
    g.copy_from(rhs);
    std::cerr<<1<<std::endl;
    return g;
  }

}

#endif

