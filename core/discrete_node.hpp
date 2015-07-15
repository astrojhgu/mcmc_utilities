#ifndef DISCRETE_NODE
#define DISCRETE_NODE
#include "stochastic_node.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class discrete_node
    :public stochastic_node<T_p,T_var1>
  {
  public:
    discrete_node(int nparents,T_var1 v_)
      :stochastic_node<T_p,T_var1>(nparents,v_)
    {}
    
    discrete_node()=delete;
    discrete_node(const discrete_node<T_p,T_var1>& )=delete;
    discrete_node<T_p,T_var1>& operator=(const discrete_node<T_p,T_var1>&)=delete;
    
  public:
    void do_sample(const base_urand<T_p>& rnd)override
    {
      this->set_value(discrete_sample(*this,rnd));
    }
  };
}



#endif
