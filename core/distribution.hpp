#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP
#include <cmath>
#include <limits>
#include "mcmc_traits.hpp"
#include "mcmc_exception.hpp"

namespace mcmc_utilities
{
  int dbg_status=0;
  template <typename T_p,typename T_var>
  class probability_density_md
  {
  public:
    int status;
    T_p eval_log(const T_var& x)const
    {
      return do_eval_log(x);
    }

    probability_density_md* clone()const
    {
      return do_clone();
    }
    
    void destroy()
    {
      do_destroy();
    }

    virtual ~probability_density_md()
    {}
    
    void var_range(T_var& x1,T_var& x2)const
    {
      do_var_range(x1,x2);
    }

  private:
    virtual T_p do_eval_log(const T_var& x)const=0;
    virtual probability_density_md* do_clone()const=0;
    virtual void do_var_range(T_var& x1,T_var& x2)const=0;
    virtual void do_destroy()
    {
      delete this;
    }
  };
}
#endif
