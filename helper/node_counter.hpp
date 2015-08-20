#ifndef NODE_COUNTER_HPP
#define NODE_COUNTER_HPP
#include <string>
#include <atomic>

namespace mcmc_utilities
{
  
  template <typename T>
  std::string node_count()
  {
    static std::atomic<int> c(0);
    ++c;
    return std::to_string(c);
  }

}

#endif
