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
  public:
    deterministic_node(size_t nparents,size_t ndim)
      :node<T,T_vector>(nparents,ndim)
    {}

    deterministic_node(size_t nparents)
      :node<T,T_vector>(nparents,1)
    {}
    
    deterministic_node()=delete;
    deterministic_node(const deterministic_node<T,T_vector>& )=delete;
    deterministic_node<T,T_vector>& operator=(const deterministic_node<T,T_vector>&)=delete;

  public:
    T calc(size_t idx,const T_vector<T>& parents)const
    {
      return this->do_calc(idx,parents);
    }

  private:
    void do_connect_to_parent(node<T,T_vector>*  rhs,size_t n,size_t idx) override
    {
      set_element(this->parents,n,std::make_pair(rhs,idx));
      rhs->add_deterministic_child(this);
    }

    virtual T do_calc(size_t idx,const T_vector<T>& parents)const=0;

    T do_value(size_t idx)const override
    {
#ifndef USE_NON_RECURSIVE
      T_vector<T> p(get_size(this->parents));
      for(size_t i=0;i<get_size(p);++i)
	{
	  set_element(p,i,this->parent(i));
	}
      return do_calc(idx,p);
#else
      std::stack<std::pair<const deterministic_node<T,T_vector>*,size_t> > node_stack;
      std::stack<size_t> leaf_num_stack;
      std::stack<T> operand_stack;
      
      node_stack.push(std::pair<const deterministic_node<T,T_vector>*,size_t>(this,idx));
      leaf_num_stack.push(0);
      for(;;)
	{
	  size_t nparents=get_size(node_stack.top().first->parents);
	  if(leaf_num_stack.top()==nparents)
	    {
	      T_vector<T> p(get_size(node_stack.top().first->parents));
	      
	      for(auto i=std::rbegin(p);i!=std::rend(p);++i)
		{
		  (*i)=operand_stack.top();
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
	      auto ptr_leaf=get_element(node_stack.top().first->parents,leaf_num_stack.top()).first;
	      size_t n=get_element(node_stack.top().first->parents,leaf_num_stack.top()).second;
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
#endif
    }
  };
}


#endif
