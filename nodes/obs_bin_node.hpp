#ifndef OBS_BIN_NODE_HPP
#define OBS_BIN_NODE_HPP
#include <core/observed_node1d.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <math/distributions.hpp>
#include <math/functions.hpp>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class obs_bin_node
    :public observed_node1d<T_p,T_var1>
  {
  public:
    obs_bin_node(const std::vector<T_var1>& d)
      :observed_node1d<T_p,T_var1>(2,d)
    {}

  public:
    T_p do_log_prior_prob()const override
    {
      T_p result=0;
      for(size_t i=0;i<this->nobs();++i)
	{
	  result+=logdbin(this->value(0,i),this->parent(0,i),this->parent(1,i));
	}
      return result;
    }
  };

  template <typename T_p,typename T_var1>
  class obs_bin_node_factory
    :public abstract_node_factory<T_p,T_var1>
  {
  public:
    obs_bin_node_factory()
      :abstract_node_factory<T_p,T_var1>({"p","n"},{"x"},{},{})
    {}
    
  public:
    std::shared_ptr<node<T_p,T_var1> >
    do_get_node(
	     const std::vector<T_var1>& scalar_param,
	     const std::vector<std::vector<T_var1> >& vector_param)const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new obs_bin_node<T_p,T_var1>(vector_param[0]));
    }      
  };

  template <typename T_p,typename T_var1>
  class obs_bin_vnode
    :public vnode<T_p,T_var1>
  {
  private:
    std::vector<T_var1> data;
  public:
    
    obs_bin_vnode(std::string n,const std::vector<T_var1>& d,
		  const std::pair<const vnode<T_p,T_var1>&,size_t>& p1,
		  const std::pair<const vnode<T_p,T_var1>&,size_t>& p2)
      :vnode<T_p,T_var1>("obs_bin",n,{p1,p2}),data(d)
    {
      this->binded=true;
    }

    obs_bin_vnode(const std::vector<T_var1>& d,
		  const std::pair<const vnode<T_p,T_var1>&,size_t>& p1,
		  const std::pair<const vnode<T_p,T_var1>&,size_t>& p2)
      :vnode<T_p,T_var1>("obs_bin",std::string("obs_bin")+node_count<obs_bin_vnode<T_p,T_var1> >(),{p1,p2}),data(d)
    {
      this->binded=true;
      this->named=false;
    }
    
    
    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new obs_bin_node<T_p,T_var1>(data));
    }
    
    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new obs_bin_vnode<T_p,T_var1>(*this));
    }
  };
  
  using obs_vbin=obs_bin_vnode<double,double>;
}


#endif
