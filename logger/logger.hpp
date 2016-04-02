#ifndef LOGGER_HPP
#define LOGGER_HPP
#include <vector>
#include <atomic>

namespace mcmc_utilities
{
  extern struct _logger
  {    
    std::atomic<int> counter;
    _logger();

    void incr();
    void decr();

    int get()const;    
  }logger;
}

#endif
