#include <core/distribution.hpp>
#include <utility>
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

    std::pair<double,double> do_var_range(const T_var& x,size_t ndim)const
    {
      return get_limits(x,ndim);
    }

    virtual pair<double,double> get_limits(const T_var& x,size_t ndim)const
    {
      return make_pair(0.,0.);
    }

    void display_limits(const T_var& x,size_t ndim)const
    {
      double x1(0),x2(0);
      std::pair<double,double> range=this->do_var_range(x,ndim);
      x1=range.first;
      x2=range.second;
      std::cout<<x1<<" "<<x2<<std::endl;
    }

    ~dist_md(){}
  };
};
