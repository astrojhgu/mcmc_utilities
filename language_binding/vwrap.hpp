#include <core/distribution.hpp>
#include <iostream>
namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  class dist_md
    :public probability_density_md<T_p,T_var>
  {
  public:
    //typedef typename element_type_trait<T_var>::element_type T_var1;
    dist_md(){}
  public:
    virtual T_p do_eval_log(const T_var& x)const=0;

    void do_var_range(double& x1,double& x2,const T_var& x,size_t ndim)const
    {
      //x1=5;
      //x2=10;
      x1=get_lower_limit(x,ndim);
      x2=get_upper_limit(x,ndim);
    }

    virtual double get_lower_limit(const T_var& x,size_t ndim)const
    {
      return 0;
    }

    virtual double get_upper_limit(const T_var& x,size_t ndim)const
    {
      return 0;
    }

    void display_limits(const T_var& x,size_t ndim)const
    {
      double x1(0),x2(0);
      this->do_var_range(x1,x2,x,ndim);
      std::cout<<x1<<" "<<x2<<std::endl;
    }

    ~dist_md(){}
  };
};
