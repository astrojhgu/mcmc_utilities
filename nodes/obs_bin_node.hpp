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
      for(int i=0;i<this->nobs();++i)
	{
	  result+=logdbin(this->value(0,i),this->parent(0,i),this->parent(1,i));
	}
      return result;
    }
  };

  template <typename T_p,typename T_var1>
  class _obs_bin_vnode
    :public _vnode<T_p,T_var1>
  {
  private:
    std::vector<T_var1> data;
  public:
    _obs_bin_vnode(std::string n,const std::vector<T_var1>& d,const std::initializer_list<std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> >& p)
      :_vnode<T_p,T_var1>("obs_bin",n,p),data(d)
    {
      this->binded=true;
    }

    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new obs_bin_node<T_p,T_var1>(data));
    }
  };

  template <typename T_p,typename T_var1>
  auto obs_bin_vnode(std::string n,const std::vector<T_var1>& d,const std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> & p1,const std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> & p2)
  {
    return shared_ptr<_vnode<T_p,T_var1> >(new _obs_bin_vnode<T_p,T_var1>(n,d,{p1,p2}));
  }

  template <typename T_p,typename T_var1>
  auto obs_bin_vnode(const std::vector<T_var1>& d,const std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> & p1,const std::pair<std::shared_ptr<_vnode<T_p,T_var1> >,size_t> & p2)
  {
    return shared_ptr<_vnode<T_p,T_var1> >(new _obs_bin_vnode<T_p,T_var1>(std::string("obs_bin")+node_count<_obs_bin_vnode<T_p,T_var1> >(),d,{p1,p2}));
  }

}


#endif
