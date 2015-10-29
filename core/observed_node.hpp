#ifndef OBSERVED_NODE
#define OBSERVED_NODE
#include <memory>
#include <vector>
#include <map>
#include "mcmc_exception.hpp"
#include "discrete_sample.hpp"
#include "base_urand.hpp"
#include "node.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class observed_node
    :public node<T_p,T_var1>
  {
  public:
    observed_node(size_t nparents,size_t ndim)
      :node<T_p,T_var1>(nparents,ndim)
    {}

    observed_node(size_t nparents)
      :node<T_p,T_var1>(nparents,1)
    {}

    observed_node()=delete;
    observed_node(const observed_node<T_p,T_var1>& )=delete;
    observed_node<T_p,T_var1>& operator=(const observed_node<T_p,T_var1>&)=delete;

  public:
    T_p log_posterior_prob()const
    {
      return log_prob()+this->log_likelihood();
    }

    T_p log_prob()const
    {
      return do_log_prob();
    }
    
    size_t nobs()const
    {
      return do_nobs();
    }



  private:
    virtual T_p do_log_prob()const=0;
    //virtual T_p do_log_prob_single(size_t obsid)const=0;

    virtual size_t do_nobs()const=0;
    
    void do_connect_to_parent(node<T_p,T_var1>*  rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
      rhs->add_observed_child(this);
    }
  };  
}


#endif
