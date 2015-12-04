#ifndef DETERMINISTIC_NODE_HPP
#define DETERMINISTIC_NODE_HPP
#include "node.hpp"

namespace mcmc_utilities
{
  template <typename T>
  class deterministic_node
    :public node<T>
  {
  public:
    deterministic_node(size_t nparents,size_t ndim)
      :node<T>(nparents,ndim)
    {}

    deterministic_node(size_t nparents)
      :node<T>(nparents,1)
    {}
    
    deterministic_node()=delete;
    deterministic_node(const deterministic_node<T>& )=delete;
    deterministic_node<T>& operator=(const deterministic_node<T>&)=delete;

  private:
    void do_connect_to_parent(node<T>*  rhs,size_t n,size_t idx) override
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
      rhs->add_deterministic_child(this);
    }
  };
}


#endif
