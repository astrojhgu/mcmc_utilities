#ifndef DETERMINISTIC_NODE_HPP
#define DETERMINISTIC_NODE_HPP
#include "node.hpp"
#include <stack>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class deterministic_node
    :public node<T,T_vector>
  {
  private:
    bool _use_system_stack;
  public:
    deterministic_node(size_t nparents,size_t ndim)
      :node<T,T_vector>(nparents,ndim),_use_system_stack(true)
    {}

    deterministic_node(size_t nparents)
      :node<T,T_vector>(nparents,1),_use_system_stack(true)
    {}
    
    deterministic_node()=delete;
    deterministic_node(const deterministic_node<T,T_vector>& )=delete;
    deterministic_node<T,T_vector>& operator=(const deterministic_node<T,T_vector>&)=delete;

  public:
    T calc(size_t idx,const T_vector<T>& parents)const
    {
      return this->do_calc(idx,parents);
    }

    void use_system_stack(bool s)
    {
      _use_system_stack=s;
    }

    bool use_system_stack()const
    {
      return _use_system_stack;
    }

  private:
    void do_connect_to_parent(node<T,T_vector>*  rhs,size_t n,size_t idx) override
    {
      this->set_parent(n,std::make_pair(rhs,idx));
      rhs->add_deterministic_child(this);
    }

    virtual T do_calc(size_t idx,const T_vector<T>& parents)const=0;
    
    T do_value(size_t idx)const override
    {
      if(use_system_stack())
	{
	  /*
	    T_vector<T> p(this->num_of_parents());
	    for(size_t i=0;i<get_size(p);++i)
	    {
	    set_element(p,i,this->parent(i));
	    }
	  */
	  auto p=this->parent_values();
	  return do_calc(idx,p);
	}
      else
	{
	  std::stack<std::pair<const deterministic_node<T,T_vector>*,size_t> > node_stack;
	  std::stack<size_t> leaf_num_stack;
	  std::stack<T> operand_stack;
	  
	  node_stack.push(std::pair<const deterministic_node<T,T_vector>*,size_t>(this,idx));
	  leaf_num_stack.push(0);
	  for(;;)
	    {
	      size_t nparents=node_stack.top().first->num_of_parents();
	      if(leaf_num_stack.top()==nparents)
		{
		  T_vector<T> p(node_stack.top().first->num_of_parents());
		  /*
		  for(auto i=std::rbegin(p);i!=std::rend(p);++i)
		    {
		      (*i)=operand_stack.top();
		      operand_stack.pop();
		    }
		  */
		  //for_each(std::rbegin(p),std::rend(p),[&operand_stack](auto& x){x=operand_stack.top();operand_stack.pop();});
		  auto size_of_p=get_size(p);
		  for(size_t i=0;i<size_of_p;++i)
		    {
		      set_element(p,size_of_p-i-1,operand_stack.top());
		      operand_stack.pop();
		    }
		  
		  operand_stack.push(node_stack.top().first->do_calc(node_stack.top().second, p));
		  node_stack.pop();
		  leaf_num_stack.pop();
		  
		  if(node_stack.empty())
		    {
		      return operand_stack.top();
		    }
		  continue;
		}
	      else
		{
		  //auto ptr_leaf=get_element(node_stack.top().first->parents,leaf_num_stack.top()).first;
		  auto ptr_leaf=(node_stack.top().first->get_parent(leaf_num_stack.top())).first;
		  size_t n=(node_stack.top().first->get_parent(leaf_num_stack.top())).second;
		  auto ptr_det_leaf=dynamic_cast<const deterministic_node<T,T_vector>*>(ptr_leaf);
		  ++leaf_num_stack.top();
		  if(ptr_det_leaf!=nullptr)
		    {
		      node_stack.push(std::pair<const deterministic_node<T,T_vector>*,size_t>(ptr_det_leaf,n));
		      leaf_num_stack.push(0);
		      continue;
		    }
		  else
		    {
		      operand_stack.push(ptr_leaf->value(n));
		      continue;
		    }
		}
	    }
	}
    }
  };
}

#endif
