#ifndef MIXTURE_NODE_HPP
#define MIXTURE_NODE_HPP
#include "stochastic_node.hpp"
#include <algorithm>
#include "mcmc_traits.hpp"
#include <memory>
#include <map>


namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class mixture_node
    :public stochastic_node<T,T_vector>
  {
  public:
    T_vector<std::shared_ptr<stochastic_node<T,T_vector> > > sub_models;
    T_vector<size_t> input_map_model;
    T_vector<size_t> parent_offsets;
  public:
    mixture_node(const T_vector<std::shared_ptr<stochastic_node<T,T_vector> > >& sm)
      :stochastic_node<T,T_vector>(std::accumulate(sm.begin(),sm.end(),size_t(0),[](size_t x,const std::shared_ptr<stochastic_node<T,T_vector> >& p){return x+p->num_of_parents();})+1,sm.front()->num_of_dims()),
      sub_models(sm),input_map_model(),parent_offsets()
    {
      for(size_t i=0;i<sub_models.size();++i)
	{
	  if(i==0)
	    {
	      parent_offsets.push_back(0);
	    }
	  else
	    {
	      parent_offsets.push_back(parent_offsets.back()+sub_models[i-1]->num_of_parents());
	    }
	  for(size_t j=0;j<sub_models[i]->num_of_parents();++j)
	    {
	      input_map_model.push_back(i);
	    }
	  
	}
    }

    T_vector<std::shared_ptr<stochastic_node<T,T_vector> > > convert_ptr(const T_vector<stochastic_node<T,T_vector>* >& p)
    {
      T_vector<std::shared_ptr<stochastic_node<T,T_vector> > > result;
      for(auto& i:p)
	{
	  result.push_back(std::shared_ptr<stochastic_node<T,T_vector> >(i));
	}
      return result;
    }
    
    //mixture_node(const T_vector<stochastic_node<T,T_vector>* >& sm)
    //:mixture_node(convert_ptr(sm))
    //{}
    
  public:
    void do_connect_to_parent(node<T,T_vector>* rhs,size_t n,size_t idx) override
    {
      this->set_parent(n,std::make_pair(rhs,idx));
      rhs->add_stochastic_child(this);
      if(n>0)
	{
	  this->sub_models.at(input_map_model.at(n-1))->set_parent(n-1-parent_offsets.at(input_map_model.at(n-1)),std::make_pair(rhs,idx));
	}
    }

    T do_log_prob()const override
    {
      size_t n=this->parent(0);
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  const_cast<stochastic_node<T,T_vector>&>(*sub_models.at(n)).set_value(i,this->value(i));
	}
      return sub_models.at(n)->log_prob();
    }

    bool is_continuous(size_t i)const override
    {
      throw mcmc_exception("should never be called");
      return true;
    }

    std::pair<T,T> do_var_range()const override
    {
      throw mcmc_exception("should never be called");
      size_t n=this->parent(0);
      return std::pair<T,T>();
    }

    void do_sample(base_urand<T>& urand)override
    {
      size_t n=this->parent(0);
      for(size_t i=0;i<sub_models.at(n)->num_of_dims();++i)
	{
	  sub_models.at(n)->set_observed(i,this->is_observed(i));
	}
      
      sub_models.at(n)->sample(urand);
      for(size_t i=0;i<sub_models.at(n)->num_of_dims();++i)
	{
	  this->set_value(i,sub_models.at(n)->value(i));
	}
    }

    void do_init_value(size_t n)override
    {
      for(size_t i=0;i<sub_models.size();++i)
	{
	  sub_models.at(n)->init_value(n);
	}
    }

    std::shared_ptr<node<T,T_vector> > do_clone()const override
    {
      T_vector<std::shared_ptr<stochastic_node<T,T_vector> > > sub_models1;
      for(auto& i:sub_models)
	{
	  sub_models1.push_back(std::dynamic_pointer_cast<stochastic_node<T,T_vector> >(i->clone()));
	}
      auto p=new mixture_node<T,T_vector>(sub_models1);
      for(size_t i=0;i<this->num_of_dims();++i)
	{
	  p->set_observed(i,this->is_observed(i));
	  p->set_value(i,this->value(i));
	}
      return std::shared_ptr<node<T,T_vector> >(p);
    }
  };
}

#endif
