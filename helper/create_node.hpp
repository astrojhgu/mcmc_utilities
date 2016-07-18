#include <memory>
#include <core/node.hpp>

namespace mcmc_utilities
{

  template <typename Tnode,template <typename> class T_vector,typename ...Args>
  std::shared_ptr<node<double,T_vector> > create_node(Args ...args)
  {
    return std::shared_ptr<node<double,T_vector> >(new Tnode(args...));
  }

}
