#include "logger.hpp"
#include <signal.h>
#include <iostream>
using namespace std;
using namespace mcmc_utilities;

namespace mcmc_utilities
{
  static void sig_handler(int sig)
  {
    std::cerr<<logger.counter<<std::endl;
    exit(0);
  }
  
  _logger::_logger()
    :counter(0)
  {
    for(int i=1;i<20;++i)
      {
	signal(i,sig_handler);
      }
  }
  
  void _logger::incr()
  {
    ++counter;
  }
  
  void _logger::decr()
  {
    --counter;
  }
  
  int _logger::get()const
  {
    return counter;
  }
  
  _logger logger;
}
