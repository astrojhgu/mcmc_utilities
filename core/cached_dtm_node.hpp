#ifndef CACHED_DTM_NODE_HPP
#define CACHED_DTM_NODE_HPP
#include "deterministic_node.hpp"
#include <stack>

namespace mcmc_utilities
{
  template <typename T>
  class cached_dtm_node
    :public deterministic_node<T>
  {
  private:
    std::vector<T> cached_parents;
    std::vector<T> cached_value;
    
  public:
    cached_dtm_node(size_t nparents,size_t ndim)
      :deterministic_node<T>(nparents,ndim),cached_parents(nparents),cached_value(ndim)
    {}
    
    cached_dtm_node(size_t nparents)
      :deterministic_node<T>(nparents,1),cached_parents(nparents),cached_value(1)
    {}
    
    cached_dtm_node()=delete;
    cached_dtm_node(const cached_dtm_node<T>& )=delete;
    cached_dtm_node<T>& operator=(const cached_dtm_node<T>&)=delete;

    T do_value(size_t idx)const override
    {
      std::vector<T> p(this->parents.size());
      bool parents_changed=false||(p.size()==0);
      for(size_t i=0;i<p.size();++i)
	{
	  p[i]=this->parent(i);
	  if(p[i]!=cached_parents.at(i))
	    {
	      parents_changed=true;
	      const_cast<cached_dtm_node<T>*>(this)->cached_parents.at(i)=p[i];
	    }
	  
	}
      
      if(parents_changed)
	{
	  double y=this->calc(idx,p);
	  const_cast<cached_dtm_node<T>*>(this)->cached_value.at(idx)=y;
	  return y;
	}
      else
	{
	  return cached_value.at(idx);
	}
    }
  };
}


#endif
