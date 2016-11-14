#include <nodes/const_node.hpp>

#ifndef DIFFERENTIABLE_DTM_NODE_HPP
#define DIFFERENTIABLE_DTM_NODE_HPP
#include "cached_dtm_node.hpp"
#include <set>
#include <stack>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class const_node;
  
  template <typename T,template <typename TE> class T_vector>
  class differentiable_dtm_node
    :public cached_dtm_node<T,T_vector>
  {
  public:
    differentiable_dtm_node(size_t nparents,size_t ndim)
      :cached_dtm_node<T,T_vector>(nparents,ndim)
    {}
    
    differentiable_dtm_node(size_t nparents)
      :cached_dtm_node<T,T_vector>(nparents,1)
    {}
    
    differentiable_dtm_node()=delete;
    differentiable_dtm_node(const differentiable_dtm_node<T,T_vector>& )=delete;
    differentiable_dtm_node<T,T_vector>& operator=(const differentiable_dtm_node<T,T_vector>&)=delete;
  protected:
    std::set<std::pair<stochastic_node<T,T_vector>*,size_t> > enumerate_stochastic_parents()const
    {
      std::set<std::pair<stochastic_node<T,T_vector>*,size_t> > result;
      for(size_t i=0;i<this->num_of_parents;++i)
      //for(auto& p:this->parents)
	{
	  auto& p=this->get_parent(i).first;
	  auto ps=dynamic_cast<stochastic_node<T,T_vector>*>(p.first);
	  auto pd=dynamic_cast<deterministic_node<T,T_vector>*>(p.first);
	  
	  if(ps==nullptr&&pd!=nullptr)
	    {
	      auto s=pd->enumerate_stochastic_parents();
	      result.insert(std::begin(s),std::end(s));
	    }
	  else if(ps!=nullptr&&pd==nullptr)
	    {
	      result.insert(p);
	    }
	  else
	    {
	      throw mcmc_exception("neither stochastic nor deterministic node or both stochastic and deterministic node");
	    }
	}
      return result;
    }

    
  };
}


#endif
