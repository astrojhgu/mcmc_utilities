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
    deterministic_node(int nparents)
      :node<T_p,T_var1>(nparents)
    {}
    
    deterministic_node()=delete;
    deterministic_node(const deterministic_node<T_p,T_var1>& )=delete;
    deterministic_node<T_p,T_var1>& operator=(const deterministic_node<T_p,T_var1>&)=delete;

  private:
    void do_connect_to_parent(node<T_p,T_var1>*  rhs,int n) override
    {
      this->parents.at(n)=rhs;
      rhs->deterministic_children.push_back(this);
    }
  };
}


#endif
