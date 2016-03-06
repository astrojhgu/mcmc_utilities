#include "error_handler.hpp"
#include "stochastic_node.hpp"

#ifndef NODE_HPP
#define NODE_HPP
#include <stack>
#include <memory>
#include <vector>
#include <list>
#include <map>

namespace mcmc_utilities
{
  template <typename T>
  class stochastic_node;

  template <typename T>
  class deterministic_node;
  
  template <typename T>
  class node
  {
    //protected:
  public:
    std::list<stochastic_node<T>* > stochastic_children;
    std::list<deterministic_node<T>* > deterministic_children;
    std::vector<std::pair<node<T>*,size_t> > parents;
    size_t ndims;
    std::vector<int> initialized;
    
  public:
    node(size_t nparents,size_t ndim1)
      :parents(nparents),ndims(ndim1),initialized(ndim1)
    {}
    
    node()=delete;
    node(const node<T>& rhs)=delete;
    node<T>& operator=(const node<T>& rhs)=delete;
    
    virtual ~node(){}

    std::shared_ptr<node<T> > clone()const
    {
      return do_clone();
    }
    
  public:
    bool is_initialized(size_t n)const
    {
      return initialized[n]!=0;
    }

    void set_initialized(size_t n,bool i)
    {
      initialized[n]=i;
    }

    void initialize()
    {
      for(size_t i=0;i<ndims;++i)
	{
	  if(!is_initialized(i))
	    {
	      initialize(i);
	    }
	}
    }

    void initialize(size_t n)
    {
      if(n<ndims)
	{
	  for(auto& p:parents)
	    {
	      for(size_t i=0;i<p.first->num_of_dims();++i)
		{
		  if(!p.first->is_initialized(i))
		    {
		      p.first->initialize(i);
		    }
		}
	    }
	  do_initialize(n);
	  initialized[n]=1;
	}
    }

    size_t num_of_parents()const
    {
      return parents.size();
    }

    const std::pair<node<T>*,size_t>& get_parent(size_t i)const
    {
      return parents[i];
    }

    size_t num_of_dims()const
    {
      return ndims;
    }
    
    T value(size_t idx)const
    {
      return this->do_value(idx);
    }

    T log_likelihood()const
    {
#ifndef USE_NON_RECURSIVE
      T result=0;
      for(auto& p : stochastic_children)
	{
	  result+=p->log_prob();
	}
      for(auto& p:deterministic_children)
	{
	  result+=p->log_likelihood();
	}
      return result;
#else
      T result=0;
      for(auto& p : stochastic_children)
	{
	  result+=p->log_prob();
	}

      std::stack<const deterministic_node<T>*> node_stack;
      std::stack<typename std::list<deterministic_node<T>* >::const_iterator> current_node_stack;
      
      for(auto& p:deterministic_children)
	{
	  node_stack.push(p);
	  current_node_stack.push(node_stack.top()->deterministic_children.begin());
	  for(;;)
	    {
	      if(current_node_stack.top()==node_stack.top()->deterministic_children.end())
		{
		  for(auto& p : node_stack.top()->stochastic_children)
		    {
		      result+=p->log_prob();
		    }
		  node_stack.pop();
		  current_node_stack.pop();
		  if(node_stack.empty())
		    {
		      break;
		    }
		}
	      else
		{
		  node_stack.push(*(current_node_stack.top()++));
		  current_node_stack.push(node_stack.top()->deterministic_children.begin());
		}
	    }
	}
      return result;
#endif
    }

    void connect_to_parent(node<T>* prhs,size_t n,size_t idx)
    {
      if(idx>=prhs->num_of_dims())
	{
	  throw output_num_mismatch();
	}
      if(n>=parents.size())
	{
	  throw parent_num_mismatch();
	}
      do_connect_to_parent(prhs,n,idx);
    }

    void add_stochastic_child(stochastic_node<T>* prhs)
    {
      stochastic_children.push_back(prhs);
    }

    void add_deterministic_child(deterministic_node<T>* prhs)
    {
      deterministic_children.push_back(prhs);
    }

    T parent(size_t pid)const
    //pid:parent id
    //obsid:the id in a set of observed values
    {
      return parents[pid].first->value(parents[pid].second);
    }
  private:
    virtual T do_value(size_t idx)const=0;
    
    virtual void do_connect_to_parent(node<T>* prhs,size_t n,size_t idx)=0;

    virtual void do_initialize(size_t n)
    {}

    virtual std::shared_ptr<node<T> > do_clone()const
    {
      not_implemented e;
      e.attach_message("clone operation not implemented");
      throw e;
      return std::shared_ptr<node<T> >(nullptr);
    }
  };
}

#endif

