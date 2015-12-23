#ifndef DETERMINISTIC_NODE_HPP
#define DETERMINISTIC_NODE_HPP
#include "node.hpp"
#include <stack>

namespace mcmc_utilities
{
  template <typename T>
  class deterministic_node
    :public node<T>
  {
  public:
    deterministic_node(size_t nparents,size_t ndim)
      :node<T>(nparents,ndim)
    {}

    deterministic_node(size_t nparents)
      :node<T>(nparents,1)
    {}
    
    deterministic_node()=delete;
    deterministic_node(const deterministic_node<T>& )=delete;
    deterministic_node<T>& operator=(const deterministic_node<T>&)=delete;


  private:
    void do_connect_to_parent(node<T>*  rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
      rhs->add_deterministic_child(this);
    }

    virtual T do_calc(size_t idx,const std::vector<T>& parents)const=0;

    T do_value(size_t idx)const override
    {
#ifdef RECURSIVE
      std::vector<T> p(this->parents.size());
      for(int i=0;i<p.size();++i)
	{
	  p[i]=this->parent(i);
	}
      return do_calc(idx,p);
#else
      std::stack<std::pair<const deterministic_node<T>*,size_t> > node_stack;
      std::stack<size_t> leaf_num_stack;
      std::stack<T> operand_stack;
      
      node_stack.push(std::pair<const deterministic_node<T>*,size_t>(this,idx));
      leaf_num_stack.push(0);
      for(;;)
	{
	  int nparents=node_stack.top().first->parents.size();
	  if(leaf_num_stack.top()==nparents)
	    {
	      std::vector<T> p(node_stack.top().first->parents.size());
	      
	      for(auto i=p.rbegin();i!=p.rend();++i)
		{
		  (*i)=operand_stack.top();
		  operand_stack.pop();
		}

	      operand_stack.push(
				      node_stack.top().first->calc(node_stack.top().second,
								   p)
				      );
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
	      auto ptr_leaf=node_stack.top().first->parents.at(leaf_num_stack.top()).first;
	      size_t n=node_stack.top().first->parents[leaf_num_stack.top()].second;
	      auto ptr_det_leaf=dynamic_cast<const deterministic_node<T>*>(ptr_leaf);
	      if(ptr_det_leaf!=nullptr)
		{
		  ++leaf_num_stack.top();
		  node_stack.push(std::pair<const deterministic_node<T>*,size_t>(ptr_det_leaf,n));
		  leaf_num_stack.push(0);
		  continue;
		}
	      else
		{
		  operand_stack.push(ptr_leaf->value(n));
		  ++leaf_num_stack.top();
		  continue;
		}
	    }
	}
#endif
    }
  };
}


#endif
