#ifndef DETERMINISTIC_NODE_HPP
#define DETERMINISTIC_NODE_HPP
#include "node.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class deterministic_node
    :public node<T_p,T_var1>
  {
  public:
    deterministic_node(size_t nparents,size_t ndim)
      :node<T_p,T_var1>(nparents,ndim)
    {}

    deterministic_node(size_t nparents)
      :node<T_p,T_var1>(nparents,1)
    {}
    
    deterministic_node()=delete;
    deterministic_node(const deterministic_node<T_p,T_var1>& )=delete;
    deterministic_node<T_p,T_var1>& operator=(const deterministic_node<T_p,T_var1>&)=delete;

  private:
    void do_connect_to_parent(node<T_p,T_var1>*  rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
      rhs->add_deterministic_child(this);
    }
  };
}


#endif
