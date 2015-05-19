#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP
#include <cmath>
#include <limits>
#include "mcmc_traits.hpp"
#include "mcmc_exception.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  class probability_density_md
  {
  public:
    typedef typename element_type_trait<T_var>::element_type T_var1;
  public:
    T_p eval_log(const T_var& x)const
    {
      return do_eval_log(x);
    }

    virtual ~probability_density_md()
    {}
    
    void var_range(T_var1& x1,T_var1& x2,const T_var& x,size_t ndim)const
    {
      do_var_range(x1,x2,x,ndim);
    }

  private:
    virtual T_p do_eval_log(const T_var& x)const=0;
    virtual void do_var_range(T_var1& x1,T_var1& x2,const T_var& x,size_t ndim)const=0;
  };
  
  template <typename T_p,typename T_var>
  class probability_density_1d
  {
  public:
    T_p eval_log(const T_var& x)const
    {
      return do_eval_log(x);
    }

    virtual ~probability_density_1d()
    {}
    
    void var_range(T_var& x1,T_var& x2)const
    {
      do_var_range(x1,x2);
    }

  private:
    virtual T_p do_eval_log(const T_var& x)const=0;
    virtual void do_var_range(T_var& x1,T_var& x2)const=0;
  };

}
#endif
