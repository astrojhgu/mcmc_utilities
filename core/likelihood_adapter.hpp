#ifndef LIKELIHOOD_ADAPTER
#define LIKELIHOOD_ADAPTER

#include <memory>
#include "stochastic_node.hpp"

namespace mcmc_utilities
{
  template <typename T>
  class likelihood_adapter
    :public stochastic_node<T>
  {
  private:
    std::shared_ptr<stochastic_node<T> > psn;

  public:
    likelihood_adapter(const std::shared_ptr<stochastic_node<T> >& _psn)
      :stochastic_node<T>(_psn->num_of_parents()+_psn->num_of_dims(),std::vector<T>()),
      psn(_psn)
    {
    }

    likelihood_adapter(stochastic_node<T>* _psn)
      :stochastic_node<T>(_psn->num_of_parents()+_psn->num_of_dims(),std::vector<T>()),
      psn(_psn)
    {
    }

  public:
    T do_log_prob()const override final
    {
      for(size_t i=psn->num_of_parents();i!=this->num_of_parents();++i)
	{
	  psn->set_value(i-psn->num_of_parents(),this->parent(i));
	}
      return psn->log_prob();
    }

    void do_connect_to_parent(node<T>* rhs,size_t n,size_t idx) override final
    {
      this->parents.at(n)=std::make_pair(rhs,idx);
      if(n<psn->num_of_parents())
	{
	  psn->connect_to_parent(rhs,n,idx);
	}
      else
	{
	  rhs->add_stochastic_child(this);
	}
    }

    bool is_continuous(size_t idx)const override final
    {
      throw mcmc_exception("should never be called!");
      return true;
    }
    
    std::pair<T,T> do_var_range()const override final
    {
      throw mcmc_exception("should never be called!");
      return std::pair<T,T>();
    } 
  };
}

#endif
