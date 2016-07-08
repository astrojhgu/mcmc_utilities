#ifndef TP_AWARE_NODE_HPP
#define TP_AWARE_NODE_HPP
#include "cached_dtm_node.hpp"
#include <set>
#include <stack>

namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class tp_aware_dtm_node
    :public cached_dtm_node<T,T_vector>
  {
  private:
    
  public:
    tp_aware_dtm_node(size_t nparents,size_t ndim)
      :cached_dtm_node<T,T_vector>(nparents,ndim)
    {}
    
    tp_aware_dtm_node(size_t nparents)
      :cached_dtm_node<T,T_vector>(nparents,1)
    {}
    
    tp_aware_dtm_node()=delete;
    tp_aware_dtm_node(const tp_aware_dtm_node<T,T_vector>& )=delete;
    tp_aware_dtm_node<T,T_vector>& operator=(const tp_aware_dtm_node<T,T_vector>&)=delete;


  protected:
    std::set<std::pair<stochastic_node<T,T_vector>*,size_t> > enumerate_stochastic_parents()const
    {
      std::set<std::pair<stochastic_node<T,T_vector>*,size_t> > result;
      for(auto& p:this->parents)
	{
	  auto ps=dynamic_cast<stochastic_node<T,T_vector>*>(p.first);
	  auto pd=dynamic_cast<deterministic_node<T,T_vector>*>(p.first);
	  
	  if(ps==nullptr&&pd!=nullptr)
	    {
	      auto s=pd->enumerate_stochastic_parents();
	      result.insert(s.begin(),s.end());
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
