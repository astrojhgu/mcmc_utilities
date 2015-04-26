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
    {};

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

  class var_range_not_encloses
    :public mcmc_exception
  {
  public:
    var_range_not_encloses()
      :mcmc_exception("var range not encloses")
    {}
  };

  class p_range_not_encloses
    :public mcmc_exception
  {
  public:
    p_range_not_encloses()
      :mcmc_exception("probability density range not encloses")
    {}
  };

  class not_enveloping
    :public mcmc_exception
  {
  public:
    not_enveloping()
      :mcmc_exception("enveloping distribution not enveloping the target distribution")
    {}
  };

  class arms_exception
    :public mcmc_exception
  {
  private:
    int err_code;
  public:
    arms_exception(int n)
      :mcmc_exception("arms exception"),err_code(n)
    {
      std::cerr<<"err code="<<n<<std::endl;
    }
  };
  
}

#endif

