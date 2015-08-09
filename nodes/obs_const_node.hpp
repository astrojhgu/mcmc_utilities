#ifndef OBS_CONST_NODE_HPP
#define OBS_CONST_NODE_HPP
#include <core/observed_node1d.hpp>
#include <helper/vnode.hpp>
#include <helper/node_counter.hpp>
#include <math/distributions.hpp>
#include <math/functions.hpp>

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class obs_const_node
    :public observed_node1d<T_p,T_var1>
  {
  public:
    obs_const_node(const std::vector<T_var1>& d)
      :observed_node1d<T_p,T_var1>(0,d)
    {}

  public:
    T_p do_log_prior_prob()const override
    {
      assert(0);//should never be called
      return 0;
    }
  };

  template <typename T_p,typename T_var1>
  class obs_const_vnode
    :public vnode<T_p,T_var1>
  {
  private:
    std::vector<T_var1> data;
  public:
    obs_const_vnode(std::string n,const std::vector<T_var1>& d)
      :vnode<T_p,T_var1>("obs_const",n,{}),data(d)
    {
      this->binded=true;
    }

    std::shared_ptr<node<T_p,T_var1> > get_node()const override
    {
      return std::shared_ptr<node<T_p,T_var1> >(new obs_const_node<T_p,T_var1>(data));
    }

    std::shared_ptr<vnode<T_p,T_var1> > clone()const override
    {
      return std::shared_ptr<vnode<T_p,T_var1> >(new obs_const_vnode<T_p,T_var1>(*this));
    }
  };

}


#endif