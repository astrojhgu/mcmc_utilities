#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP
#include <cmath>
#include <limits>
#include "mcmc_traits.hpp"
#include "mcmc_exception.hpp"
#include <vector>

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

    size_t num_init_points(const T_var& x,size_t ndim)const
    {
      return do_num_init_points(x,ndim);
    }

    T_var1 init_point(size_t n,const T_var& x,size_t ndim)const
    {
      return do_init_point(n,x,ndim);
    }
      

  private:
    virtual T_p do_eval_log(const T_var& x)const=0;
    virtual void do_var_range(T_var1& x1,T_var1& x2,const T_var& x,size_t ndim)const=0;
    virtual size_t do_num_init_points(const T_var& x,size_t ndim)const
    {
      return 0;
    }
    virtual T_var1 do_init_point(size_t n,const T_var& x,size_t ndim)const
    {
      return T_var1();
    }
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

    size_t num_init_points()const
    {
      return do_num_init_points();
    }

    T_var init_point(size_t n)const
    {
      return do_init_point(n);
    }

  private:
    virtual T_p do_eval_log(const T_var& x)const=0;
    virtual void do_var_range(T_var& x1,T_var& x2)const=0;
    virtual size_t do_num_init_points()const
    {
      return 3;
    };
    virtual T_var do_init_point(size_t n)const
    {
      if(n!=1)
	{
	  T_var xl,xr;
	  var_range(xl,xr);
	  return xl+(xr-xl)/(num_init_points()+1)*(n+1);
	}
      return find_peak(*this);
    };
  };
}
#endif
