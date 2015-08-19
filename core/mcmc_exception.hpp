#ifndef MCMC_EXCEPTION
#define MCMC_EXCEPTION

#include <string>
#include <exception>
#include <iostream>
namespace mcmc_utilities
{
  class mcmc_exception
    :public std::exception
  {
  private:
    std::string _what;
  public:
    mcmc_exception()
    {}

    ~mcmc_exception()throw()
    {}

    mcmc_exception(const std::string& str)
      :_what(str)
    {}

    const char* what()const throw()
    {
      return _what.c_str();
    }
  };

  class index_out_of_range
    :public mcmc_exception
  {
  public:
    index_out_of_range()
      :mcmc_exception("index out of range")
    {}
  };



  class var_out_of_range
    :public mcmc_exception
  {
  public:
    var_out_of_range()
      :mcmc_exception("var out of range")
    {}
  };

  class no_candidate_points
    :public mcmc_exception
  {
  public:
    no_candidate_points()
      :mcmc_exception("no candidate points")
    {}
  };
  
  
}

#endif

