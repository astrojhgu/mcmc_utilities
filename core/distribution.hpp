#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP
#include <cmath>
#include <limits>
#include "mcmc_traits.hpp"
#include "error_handler.hpp"
//#include "find_peak.hpp"
#include <vector>
#include <utility>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class probability_density_1d;

  //template<typename T_p,typename T_var>
  //T_var find_peak(const probability_density_1d<T_p,T_var>& dist);

  
  template <typename T_p,typename T_var,template <typename TE> class T_vector>
  class probability_density_md
  {
  public:
    using T_var1=typename element_type_trait<T_var>::element_type;
  public:
    T_p eval_log(const T_var& x,int n=-1)const
    {
      return do_eval_log(x,n);
    }

    virtual ~probability_density_md()
    {}
    
    std::pair<T_var1,T_var1> var_range(const T_var& x,size_t ndim)const
    {
      return do_var_range(x,ndim);
    }

    T_vector<T_var1> init_points(const T_var& x,size_t ndim)const
    {
      return do_init_points(x,ndim);
    }
    
    T_vector<T_var1> candidate_points(const T_var& x,size_t ndim)const
    {
      return do_candidate_points(x,ndim);
    }

  private:
    virtual T_p do_eval_log(const T_var& x,int n)const=0;
    virtual std::pair<T_var1,T_var1> do_var_range(const T_var& x,size_t ndim)const=0;
    virtual T_vector<T_var1> do_init_points(const T_var& x,size_t ndim)const
    {
      return T_vector<T_var1>();
    }

    virtual T_vector<T_var1> do_candidate_points(const T_var& x,size_t ndim)const
    {
      return T_vector<T_var1>();
    }
    
  };
  
  template <typename T,template <typename TE> class T_vector>
  class probability_density_1d
  {
  public:
    T eval_log(const T& x)const
    {
      T result=do_eval_log(x);
      return result;
    }

    virtual ~probability_density_1d()
    {}
    
    std::pair<T,T> var_range()const
    {
      return do_var_range();
    }

    T_vector<T> init_points()const
    {
      return do_init_points();
    }

    T_vector<T> candidate_points()const
    {
      return do_candidate_points();
    }
    
  private:
    virtual T do_eval_log(const T& x)const=0;
    virtual std::pair<T,T> do_var_range()const=0;
    virtual T_vector<T> do_init_points()const
    {
      T_vector<T> result(5);
      for(size_t n=0;n<get_size(result);++n)
	{
	  //if(n!=1)
	    {
	      std::pair<T,T> xrange(var_range());
	      T xl=xrange.first,xr=xrange.second;
	      
	      set_element(result,n,xl+(xr-xl)/(get_size(result)+1)*(n+1));
	    }
	}
      return result;
    };

    //possible values when this distribution is discrete
    virtual T_vector<T> do_candidate_points()const
    {
      return T_vector<T>();
    }
  };
}
#endif
