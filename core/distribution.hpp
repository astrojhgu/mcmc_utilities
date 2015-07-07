#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP
#include <cmath>
#include <limits>
#include "mcmc_traits.hpp"
#include "mcmc_exception.hpp"
#include <vector>
#include <utility>

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
    
    std::pair<T_var1,T_var1> var_range(const T_var& x,size_t ndim)const
    {
      return do_var_range(x,ndim);
    }

    std::vector<T_var1> init_points(const T_var& x,size_t ndim)const
    {
      return do_init_points(x,ndim);
    }
      

  private:
    virtual T_p do_eval_log(const T_var& x)const=0;
    virtual std::pair<T_var1,T_var1> do_var_range(const T_var& x,size_t ndim)const=0;
    virtual std::vector<T_var1> do_init_points(const T_var& x,size_t ndim)const
    {
      return std::vector<T_var1>();
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
    
    std::pair<T_var,T_var> var_range()const
    {
      return do_var_range();
    }

    std::vector<T_var> init_points()const
    {
      return do_init_points();
    }

  private:
    virtual T_p do_eval_log(const T_var& x)const=0;
    virtual std::pair<T_var,T_var> do_var_range()const=0;
    virtual std::vector<T_var> do_init_points()const
    {
      std::vector<T_var> result(3);
      for(int n=0;n<result.size();++n)
	{
	  if(n!=1)
	    {
	      std::pair<T_var,T_var> xrange(var_range());
	      T_var xl=xrange.first,xr=xrange.second;
	      
	      result[n]= xl+(xr-xl)/(4)*(n+1);
	    }
	  else
	    {
	      result[n]= find_peak(*this);
	    }
	}
      return result;
    };
  };
}
#endif
