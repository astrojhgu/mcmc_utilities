#ifndef STOCHASTIC_NODE
#define STOCHASTIC_NODE
#include <memory>
#include <vector>
#include <map>
#include "mcmc_exception.hpp"
#include "arms.hpp"
#include "discrete_sample.hpp"
#include "base_urand.hpp"
#include "node.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class stochastic_node
    :public node<T_p,T_var1>,public probability_density_1d<T_p,T_var1>
  {
  private:
    T_p v;

  public:
    stochastic_node(int nparents,T_var1 v_)
      :node<T_p,T_var1>(nparents),
      v(v_)
    {}

    stochastic_node()=delete;
    stochastic_node(const stochastic_node<T_p,T_var1>& )=delete;
    stochastic_node<T_p,T_var1>& operator=(const stochastic_node<T_p,T_var1>&)=delete;

  public:
    T_p log_post_prob()const
    {
      return log_prior_prob()+this->log_likelihood();
    }

    T_p do_eval_log(const T_var1& x)const override
    {
      const_cast<stochastic_node*>(this)->set_value(x);
      return log_post_prob();
    }
    
    T_p log_prior_prob()const
    {
      return do_log_prior_prob();
    }

    void set_value(const T_var1& v_)
    {
      v=v_;
    }

  public:
    void sample(const base_urand<T_p>& rnd)
    {
      do_sample(rnd);
    }
  private:
    virtual T_p do_log_prior_prob()const=0;
    T_var1 do_value()const override
    {
      return v;
    }

    void do_connect_to_parent(node<T_p,T_var1>*  rhs,int n) override
    {
      this->parents.at(n)=rhs;
      rhs->stochastic_children.push_back(this);
    }

    virtual void do_sample(const base_urand<T_p>&)=0;
  };  
}


#endif
