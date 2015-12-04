#include "mcmc_exception.hpp"
#include "stochastic_node.hpp"

#ifndef NODE_HPP
#define NODE_HPP

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
  protected:
    std::list<stochastic_node<T>* > stochastic_children;
    std::list<deterministic_node<T>* > deterministic_children;
    std::vector<std::pair<node<T>*,size_t> > parents;
    size_t ndim_;
    
  public:
    node(size_t nparents,size_t ndim1)
      :parents(nparents),ndim_(ndim1)
    {}
    
    node()=delete;
    node(const node<T>& rhs)=delete;
    node<T>& operator=(const node<T>& rhs)=delete;
    
    virtual ~node(){}
    
  public:
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
      return ndim_;
    }
    
    T value(size_t idx)const
    {
      return do_value(idx);
    }

    T regulate(const T& x)const
    {
      return do_regulate(x);
    }

    virtual T log_likelihood()const final
    {
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

    virtual T do_regulate(const T& x)const
    {
      return x;
    }
  };
}

#endif

