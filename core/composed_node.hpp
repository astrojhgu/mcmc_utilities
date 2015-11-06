#ifndef COMPOSED_NODE_HPP
#define COMPOSED_NODE_HPP
#include "deterministic_node.hpp"
#include <map>


namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class forward_node
    :public deterministic_node<T_p,T_var1>
  {
  public:
    forward_node()
      :deterministic_node<T_p,T_var1>(1,1)
      {}

    T_var1 do_value(size_t idx)const override
    {
      return this->parent(0);
    }

    void do_connect_to_parent(node<T_p,T_var1>* rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
    }

  };

  
  
  template <typename T_p,typename T_var1,typename T_tag=std::string>
  class composed_node
    :public deterministic_node<T_p,T_var1>
  {
  private:
    std::map<T_tag,std::shared_ptr<deterministic_node<T_p,T_var1> > >elements;
    std::vector<std::pair<std::shared_ptr<deterministic_node<T_p,T_var1> >,size_t> > param_list;
    std::vector<std::pair<std::shared_ptr<deterministic_node<T_p,T_var1> >,size_t> > return_list;
  public:
    composed_node(size_t nparents,size_t ndim)
      :deterministic_node<T_p,T_var1>(nparents,ndim)
    {}
    
  public:
    void add_node(deterministic_node<T_p,T_var1>* pn,
		  const T_tag& tag,
		  const std::vector<std::pair<T_tag,size_t> >& parents,
		  const std::vector<size_t> rlist
		  )
    {
      add_node(shared_ptr<deterministic_node<T_p,T_var1> >(pn),
	       tag,parents,rlist);
    }
    
    void add_node(const std::shared_ptr<deterministic_node<T_p,T_var1> >& pn,
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
	      param_list.push_back({std::shared_ptr<deterministic_node<T_p,T_var1> >(new forward_node<T_p,T_var1>),i});
	      pn->connect_to_parent(param_list.back().first.get(),i,0);
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

    void do_connect_to_parent(node<T_p,T_var1>* rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
      rhs->add_deterministic_child(this);
            
      param_list[n].first->connect_to_parent(rhs,0,idx);
    }

    T_var1 do_value(size_t idx)const override
    {
      return return_list[idx].first->value(return_list[idx].second);
    }
  };
}

#endif
