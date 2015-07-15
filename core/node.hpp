#ifndef NODE_HPP
#define NODE_HPP

#include <memory>
#include <vector>
#include <map>
#include "mcmc_exception.hpp"


namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class stochastic_node;

  template <typename T_p,typename T_var1>
  class deterministic_node;
  
  template <typename T_p,typename T_var1>
  class node
  {
    friend class stochastic_node<T_p,T_var1>;
    friend class deterministic_node<T_p,T_var1>;
  protected:
    std::vector<stochastic_node<T_p,T_var1>* > stochastic_children;
    std::vector<deterministic_node<T_p,T_var1 >* > deterministic_children;
    std::vector<node<T_p,T_var1>* > parents;

  public:
    explicit node(int nparents)
      :parents(nparents)
    {}

    node()=delete;
    node(const node<T_p,T_var1>& rhs)=delete;
    node<T_p,T_var1>& operator=(const node<T_p,T_var1>& rhs)=delete;

    virtual ~node(){}
    
  public:
    int num_of_parents()const
    {
      return parents.size();
    }
    
    T_var1 value()const
    {
      return do_value();
    }

    virtual T_p log_likelihood()const final
    {
      T_p result=0;
      for(auto& p : stochastic_children)
	{
	  result+=p->log_prior_prob();
	}
      for(auto& p:deterministic_children)
	{
	  result+=p->log_likelihood();
	}
      return result;
    }

    void connect_to_parent(node<T_p,T_var1>* prhs,int n)
    {
      do_connect_to_parent(prhs,n);
    }
  private:
    virtual T_var1 do_value()const=0;
    virtual void do_connect_to_parent(node<T_p,T_var1>* prhs,int n)=0;
  };
}

#endif

