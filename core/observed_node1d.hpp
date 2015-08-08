#ifndef OBSERVED_NODE1D
#define OBSERVED_NODE1D
#include "observed_node.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class observed_node1d
    :public observed_node<T_p,T_var1>
  {
  private:
    std::vector<T_var1> data;

  public:
    observed_node1d(size_t nparents,const std::vector<T_var1>& d)
      :observed_node<T_p,T_var1>(nparents),data(d)
    {}

  private:
    size_t do_nobs()const override
    {
      return data.size();
    }

    T_var1 do_value(size_t idx,size_t i)const override
    {
      return data.at(i);
    }
  };
}


#endif
