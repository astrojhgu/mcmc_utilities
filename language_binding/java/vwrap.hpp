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

    void do_var_range(double& x1,double& x2,const T_var& x,size_t ndim)const
    {
      std::pair<double,double> result(get_limits(x,ndim));
      x1=result.first;
      x2=result.second;
    }

    virtual std::pair<double,double> get_limits(const T_var& x,size_t ndim)const
    {
      return make_pair(0.,0.);
    }

    void display_limits(const T_var& x,size_t ndim)const
    {
      double x1(0),x2(0);
      this->do_var_range(x1,x2,x,ndim);
      std::cout<<x1<<" "<<x2<<std::endl;
    }

    ~dist_md(){}
  };

  template <typename T_p,typename T_var,typename T_urand>
  void gibbs_sample_real(const dist_md<T_p,T_var>& pd,T_var& init_var,bool dometrop,const T_urand& urand, int sample_cnt=10)
  {
    //cout<<"fsfaf"<<endl;
    gibbs_sample(pd,init_var,dometrop,urand,sample_cnt);
  }

  double f()
  {
    return 10;
  }
  
};
