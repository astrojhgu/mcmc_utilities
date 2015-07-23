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
    continuous_node(size_t nparents,const std::vector<T_var1>& v_)
      :stochastic_node<T_p,T_var1>(nparents,v_)
    {}

    continuous_node(size_t nparents,T_var1 v_)
      :stochastic_node<T_p,T_var1>(nparents,v_)
    {}
    
    continuous_node()=delete;
    continuous_node(const continuous_node<T_p,T_var1>& )=delete;
    continuous_node<T_p,T_var1>& operator=(const continuous_node<T_p,T_var1>&)=delete;

  public:
    void do_sample(const base_urand<T_p>& rnd)override
    {
      constexpr size_t nsamp=10;
      for(int i=0;i<this->num_of_dims();++i)
	{
	  this->set_current_idx(i);
	  T_var1 xprev=this->value(i,0);
	  std::vector<T_var1> xsamp(nsamp);
	  arms_simple(*this,xprev,xsamp,dometrop(),rnd);
	  this->set_value(i,xsamp.back());
	}
    }
  private:
    virtual bool dometrop()const
    {
      return true;
    }
  };
}


#endif
