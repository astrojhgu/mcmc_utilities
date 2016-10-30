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
      :_what()
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

    void attach_message(const std::string& str)
    {
      _what+="\n";
      _what+=str;
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

  class no_candidate_point
    :public mcmc_exception
  {
  public:
    no_candidate_point()
      :mcmc_exception("no candidate point")
    {}
  };

  class too_few_init_points
    :public mcmc_exception
  {
  public:
    too_few_init_points()
      :mcmc_exception("too_few_init_points")
    {}
  };

  class data_not_in_order
    :public mcmc_exception
  {
  public:
    data_not_in_order()
      :mcmc_exception("data_not_in_order")
    {}
  };

  class nan_or_inf
    :public mcmc_exception
  {
  public:
    nan_or_inf()
      :mcmc_exception("nan_or_inf")
    {}    
  };

  class cum_lt_zero
    :public mcmc_exception
  {
  public:
    cum_lt_zero()
      :mcmc_exception("cum < 0")
    {}    
  };

  class more_init_points_needed
    :public mcmc_exception
  {
  public:
    more_init_points_needed()
      :mcmc_exception("more init points needed")
    {}    
  };

  class search_failed
    :public mcmc_exception
  {
  public:
    search_failed()
      :mcmc_exception("search failed")
    {}    
  };


  class node_name_already_used
    :public mcmc_exception
  {
  public:
    node_name_already_used()
      :mcmc_exception("node name already used")
    {}    
  };

  class node_already_added
    :public mcmc_exception
  {
  public:
    node_already_added()
      :mcmc_exception("node has already added")
    {}
  };

  class parent_not_connected
    :public mcmc_exception
  {
  public:
    parent_not_connected()
      :mcmc_exception("parent not connected")
    {}
  };
  
  class parents_not_exist
    :public mcmc_exception
  {
  public:
    parents_not_exist()
      :mcmc_exception("parents not exist")
    {}    
  };

  class parent_already_connected
    :public mcmc_exception
  {
  public:
    parent_already_connected()
      :mcmc_exception("parent alread connected")
    {}
  };

  class parent_num_mismatch
    :public mcmc_exception
  {
  public:
    parent_num_mismatch()
      :mcmc_exception("parent_num_mismatch")
    {}    
  };

  class output_num_mismatch
    :public mcmc_exception
  {
  public:
    output_num_mismatch()
      :mcmc_exception("output_num_mismatch")
    {}    
  };


  class node_not_found
    :public mcmc_exception
  {
  public:
    node_not_found()
      :mcmc_exception("node not found")
    {}

    node_not_found(const std::string& node_name)
      :mcmc_exception(node_name+" not found")
    {}
  };

  class invalid_node_type
    :public mcmc_exception
  {
  public:
    invalid_node_type()
      :mcmc_exception("invalid node type")
    {}    
  };

  class pointer_expired
    :public mcmc_exception
  {
  public:
    pointer_expired()
      :mcmc_exception("expired")
    {}    
  };

  class too_many_rejections
    :public mcmc_exception
  {
  public:
    too_many_rejections()
      :mcmc_exception("too may rejections, maybe illed-shaped pdf")
    {}    
  };

  class ill_conditioned_distribution
    :public mcmc_exception
  {
  public:
    ill_conditioned_distribution()
      :mcmc_exception("ill condition distribution")
    {}

    ill_conditioned_distribution(const std::string& str)
      :mcmc_exception(str)
    {}
  };

  class not_implemented
    :public mcmc_exception
  {
  public:
    not_implemented()
      :mcmc_exception("not implemented")
    {}
  };

  class too_small_var_range
    :public mcmc_exception
  {
  public:
    too_small_var_range()
      :mcmc_exception("too small var range")
    {}
  };
}

#endif

