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

  class eval_method_not_implemented
    :public mcmc_exception
  {
  public:
    eval_method_not_implemented()
      :mcmc_exception("eval method not implemented, try eval_log")
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

  class invalid_prob_value
    :public mcmc_exception
  {
  public:
    invalid_prob_value()
      :mcmc_exception("invalid probability value, probably overflow")
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

  class centile_out_of_range
    :public mcmc_exception
  {
  public:
    centile_out_of_range()
      :mcmc_exception("centile out of range")
    {}
  };

  class previous_iter_out_of_range
    :public mcmc_exception
  {
  public:
    previous_iter_out_of_range()
      :mcmc_exception("previous iteration out of range")
    {}
  };

  class envelope_error
    :public mcmc_exception
  {
  public:
    envelope_error()
      :mcmc_exception("envelope error")
    {}
  };
  
  class too_few_init_points
    :public mcmc_exception
  {
  public:
    too_few_init_points()
      :mcmc_exception("too few init point")
    {}
  };
  
  class too_many_init_points
    :public mcmc_exception
  {
  public:
    too_many_init_points()
      :mcmc_exception("too many init point")
    {}
  };

  class init_point_out_of_range
    :public mcmc_exception
  {
  public:
    init_point_out_of_range()
      :mcmc_exception("init point out of range")
    {}
  };

  class data_not_ordered
    :public mcmc_exception
  {
  public:
    data_not_ordered()
      :mcmc_exception("data not ordered")
    {}
  };
  
  class negative_convexity
    :public mcmc_exception
  {
  public:
    negative_convexity()
      :mcmc_exception("negative convexity")
    {}
  };

  class insufficient_space
    :public mcmc_exception
  {
  public:
    insufficient_space()
      :mcmc_exception("insufficient space")
    {}
  };

  class range_not_ordered
    :public mcmc_exception
  {
  public:
    range_not_ordered()
      :mcmc_exception("range not in order")
    {}
  };


  class index_out_of_range
    :public mcmc_exception
  {
  public:
    index_out_of_range()
      :mcmc_exception("index out of range")
    {}
  };

  class envelope_violation_without_metropolis
    :public mcmc_exception
  {
  public:
    envelope_violation_without_metropolis()
      :mcmc_exception("envelope violation without metropolis")
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

