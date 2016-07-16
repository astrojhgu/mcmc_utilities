#ifndef COMPOSED_NODE_HPP
#define COMPOSED_NODE_HPP
#include "cached_dtm_node.hpp"
#include "mcmc_traits.hpp"
#include <memory>
#include <map>


namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class forward_node
    :public cached_dtm_node<T,T_vector>
  {
  public:
    forward_node()
      :cached_dtm_node<T,T_vector>(1,1)
      {}

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      return get_element(parent,0);
    }

    void do_connect_to_parent(node<T,T_vector>* rhs,size_t n,size_t idx) override
    {
      //set_element(this->parents,n,std::make_pair(rhs,idx));
      this->set_parent(n,std::make_pair(rhs,idx));
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      throw mcmc_exception("should not be called!");
      return std::shared_ptr<node<T,T_vector> >(new forward_node<T,T_vector>);
    }

    T do_value(size_t idx)const override
    {
      return this->get_parent(0).first->value(idx);
    }
  };

  
  
  template <typename T,typename T_tag,template <typename TE> class T_vector>
  class composed_node
    :public cached_dtm_node<T,T_vector>
  {
  protected:
    std::map<T_tag,std::shared_ptr<deterministic_node<T,T_vector> > >elements;
    T_vector<std::shared_ptr<deterministic_node<T,T_vector> > > param_list;
    T_vector<std::pair<std::shared_ptr<deterministic_node<T,T_vector> >,size_t> > return_list;
  public:
    composed_node(size_t nparents,size_t ndim)
      :cached_dtm_node<T,T_vector>(nparents,ndim)
    {}
    
  public:
    void add_node(deterministic_node<T,T_vector>* pn,
		  const T_tag& tag,
		  //if the parent tag is not registered, then this node will registed as an input parameter
		  const T_vector<std::pair<T_tag,size_t> >& parents,
		  //which outputs whill outputed of this node
		  const T_vector<size_t> rlist
		  )
    {
      add_node(std::shared_ptr<deterministic_node<T,T_vector> >(pn),
	       tag,parents,rlist);
    }
    
    void add_node(const std::shared_ptr<deterministic_node<T,T_vector> >& pn,
		  const T_tag& tag,
		  const T_vector<std::pair<T_tag,size_t> >& parents,
		  const T_vector<size_t> rlist
		  )
    {
      if (elements.count(tag)!=0)
	{
	  throw node_name_already_used();
	}
      if(pn->num_of_parents()!=this->num_of_parents())
	{
	  throw parent_num_mismatch();
	}
      for(size_t i=0;i<this->num_of_parents();++i)
	{
	  auto iter=elements.find(this->get_parent(i).first);
	  if(iter==elements.end())
	    {
	      if(get_size(param_list)==this->num_of_parents())
		{
		  throw parent_num_mismatch();
		}
	      //param_list.push_back(make_pair(pn,i));
	      push_back(param_list,std::shared_ptr<deterministic_node<T,T_vector> >(new forward_node<T,T_vector>));
	      pn->connect_to_parent(last_element(param_list).get(),i,0);
	    }
	  else
	    {
	      pn->connect_to_parent(iter->second.get(),i,this->get_parent(i).second);
	    }
	}
      get_element(elements,tag)=pn;
      for(auto& i : rlist)
	{
	  if(get_size(return_list)==this->num_of_dims())
	    {
	      throw output_num_mismatch();
	    }
	  if(i>=pn->num_of_dims())
	    {
	      throw output_num_mismatch();
	    } 
	  push_back(return_list,std::pair<std::shared_ptr<deterministic_node<T,T_vector> >,size_t>{pn,i});
	}
    }

    void do_connect_to_parent(node<T,T_vector>* rhs,size_t n,size_t idx) override
    {
      this->set_parent(n,std::make_pair(rhs,idx));
      rhs->add_deterministic_child(this);
            
      get_element(param_list,n)->connect_to_parent(rhs,0,idx);
    }

    T do_value(size_t idx)const override
    {
      return get_element(return_list,idx).first->value(get_element(return_list,idx).second);
    }

    T do_calc(size_t idx,const T_vector<T>& parent)const override
    {
      throw mcmc_exception("should never be called!");
    }
  };
}

#endif
