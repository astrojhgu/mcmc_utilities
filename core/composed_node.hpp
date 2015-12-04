#ifndef COMPOSED_NODE_HPP
#define COMPOSED_NODE_HPP
#include "deterministic_node.hpp"
#include <memory>
#include <map>


namespace mcmc_utilities
{
  template <typename T>
  class forward_node
    :public deterministic_node<T>
  {
  public:
    forward_node()
      :deterministic_node<T>(1,1)
      {}

    T do_value(size_t idx)const override
    {
      return this->parent(0);
    }

    void do_connect_to_parent(node<T>* rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
    }

  };

  
  
  template <typename T,typename T_tag=std::string>
  class composed_node
    :public deterministic_node<T>
  {
  protected:
    std::map<T_tag,std::shared_ptr<deterministic_node<T> > >elements;
    std::vector<std::shared_ptr<deterministic_node<T> > > param_list;
    std::vector<std::pair<std::shared_ptr<deterministic_node<T> >,size_t> > return_list;
  public:
    composed_node(size_t nparents,size_t ndim)
      :deterministic_node<T>(nparents,ndim)
    {}
    
  public:
    void add_node(deterministic_node<T>* pn,
		  const T_tag& tag,
		  //if the parent tag is not registered, then this node will registed as an input parameter
		  const std::vector<std::pair<T_tag,size_t> >& parents,
		  //which outputs whill outputed of this node
		  const std::vector<size_t> rlist
		  )
    {
      add_node(std::shared_ptr<deterministic_node<T> >(pn),
	       tag,parents,rlist);
    }
    
    void add_node(const std::shared_ptr<deterministic_node<T> >& pn,
		  const T_tag& tag,
		  const std::vector<std::pair<T_tag,size_t> >& parents,
		  const std::vector<size_t> rlist
		  )
    {
      if (elements.count(tag)!=0)
	{
	  throw node_name_already_used();
	}
      if(pn->num_of_parents()!=parents.size())
	{
	  throw parent_num_mismatch();
	}
      for(size_t i=0;i<parents.size();++i)
	{
	  auto iter=elements.find(parents[i].first);
	  if(iter==elements.end())
	    {
	      if(param_list.size()==this->num_of_parents())
		{
		  throw parent_num_mismatch();
		}
	      //param_list.push_back(make_pair(pn,i));
	      param_list.push_back(std::shared_ptr<deterministic_node<T> >(new forward_node<T>));
	      pn->connect_to_parent(param_list.back().get(),i,0);
	    }
	  else
	    {
	      pn->connect_to_parent(iter->second.get(),i,parents[i].second);
	    }
	}
      elements[tag]=pn;
      for(auto& i : rlist)
	{
	  if(return_list.size()==this->num_of_dims())
	    {
	      throw output_num_mismatch();
	    }
	  if(i>=pn->num_of_dims())
	    {
	      throw output_num_mismatch();
	    } 
	  return_list.push_back({pn,i});
	}
    }

    void do_connect_to_parent(node<T>* rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
      rhs->add_deterministic_child(this);
            
      param_list[n]->connect_to_parent(rhs,0,idx);
    }

    T do_value(size_t idx)const override
    {
      return return_list[idx].first->value(return_list[idx].second);
    }
  };
}

#endif
