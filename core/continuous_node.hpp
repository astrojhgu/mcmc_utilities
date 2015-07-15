#ifndef CONTINUOUS_HPP
#define CONTINUOUS_HPP
#include "stochastic_node.hpp"

namespace mcmc_utilities
{
    template <typename T_p,typename T_var1>
  class continuous_node
    :public stochastic_node<T_p,T_var1>
  {
  public:
    continuous_node(int nparents,T_var1 v_)
      :stochastic_node<T_p,T_var1>(nparents,v_)
    {}
    
    continuous_node()=delete;
    continuous_node(const continuous_node<T_p,T_var1>& )=delete;
    continuous_node<T_p,T_var1>& operator=(const continuous_node<T_p,T_var1>&)=delete;

  public:
    void do_sample(const base_urand<T_p>& rnd)override
    {
      T_var1 xprev=this->value();
      constexpr int nsamp=10;
      std::vector<T_var1> xsamp(nsamp);
      arms_simple(*this,xprev,xsamp,dometrop(),rnd);
      this->set_value(xsamp.back());
    }
  private:
    virtual bool dometrop()const
    {
      return true;
    }
  };
}


#endif
