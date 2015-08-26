#include "mcmc_exception.hpp"
#include "stochastic_node.hpp"
#include "observed_node.hpp"

#ifndef NODE_HPP
#define NODE_HPP

#include <memory>
#include <vector>
#include <list>
#include <map>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class stochastic_node;

  template <typename T_p,typename T_var1>
  class observed_node;

  template <typename T_p,typename T_var1>
  class deterministic_node;
  
  template <typename T_p,typename T_var1>
  class node
  {
  protected:
    std::list<stochastic_node<T_p,T_var1>* > stochastic_children;
    std::list<observed_node<T_p,T_var1>* > observed_children;
    std::list<deterministic_node<T_p,T_var1 >* > deterministic_children;
    std::vector<std::pair<node<T_p,T_var1>*,size_t> > parents;
    size_t ndim_;
  public:
    node(size_t nparents,size_t ndim1)
      :parents(nparents),ndim_(ndim1)
    {}
    
    node()=delete;
    node(const node<T_p,T_var1>& rhs)=delete;
    node<T_p,T_var1>& operator=(const node<T_p,T_var1>& rhs)=delete;
    
    virtual ~node(){}
    
  public:
    size_t num_of_parents()const
    {
      return parents.size();
    }

    size_t num_of_dims()const
    {
      return ndim_;
    }
    
    T_var1 value(size_t idx,size_t obsid)const
    {
      return do_value(idx,obsid);
    }

    virtual T_p log_likelihood()const final
    {
      T_p result=0;
      for(auto& p : stochastic_children)
	{
	  result+=p->log_prior_prob();
	}
      for(auto& p: observed_children)
	{
	  result+=p->log_prior_prob();
	}
      for(auto& p:deterministic_children)
	{
	  result+=p->log_likelihood();
	}
      return result;
    }

    void connect_to_parent(node<T_p,T_var1>* prhs,size_t n,size_t idx)
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

    void add_observed_child(observed_node<T_p,T_var1>* prhs)
    {
      observed_children.push_back(prhs);
    }

    void add_stochastic_child(stochastic_node<T_p,T_var1>* prhs)
    {
      stochastic_children.push_back(prhs);
    }

    void add_deterministic_child(deterministic_node<T_p,T_var1>* prhs)
    {
      deterministic_children.push_back(prhs);
    }

    T_var1 parent(size_t pid,size_t obsid)const
    //pid:parent id
    //obsid:the id in a set of observed values
    {
      return parents[pid].first->value(parents[pid].second,obsid);
    }
  private:
    virtual T_var1 do_value(size_t idx,size_t obsid)const=0;
    
    virtual void do_connect_to_parent(node<T_p,T_var1>* prhs,size_t n,size_t idx)=0;
  };
}

#endif

