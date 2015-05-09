#ifndef DBIN
#define DBIN
#include <core/distribution.hpp>

namespace mcmc_utilities
{

  template <typename T_p,typename T_var>
  class dbin
    :public probability_density_md<T_p,T_var>
  {
  private:
    T_p p;
    T_var n;
  public:
    T_p log_factorial(const T_var& n)const
    {
      if(n==0)
	{
	  return 0;
	}
      T_p result=0;

      for(T_var x=1;x<=n;++x)
	{
	  result+=std::log(x);
	}
      return result;
    }
    T_var log_CN(const T_var& m,const T_var& n)const
    {
      //return factorial(m)/(factorial(n)*factorial(m-n));
      return log_factorial(m)-log_factorial(n)-log_factorial(m-n);
    }
  public:
    dbin(const T_p& _p,const T_var& _n)
      :p(_p),n(_n)
    {
    }
    
    dbin(const dbin& rhs)
      :p(rhs.p),n(rhs.n)
    {
    }
  private:
    T_p do_eval_log(const T_var& x)const
    {
      return log_CN(n,x)+x*std::log(p)+(n-x)*std::log(1-p);
    }
    
  dbin* do_clone()const
    {
      return new dbin(*this);
    }
    
    void do_var_range(T_var& x1,T_var& x2)const
    {
      x1=0;
      x2=n;
    }
  };
}


#endif
