#ifndef CACHED_DTM_NODE_HPP
#define CACHED_DTM_NODE_HPP
#include "deterministic_node.hpp"
#include <stack>
#include <mutex>
namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class cached_dtm_node
    :public deterministic_node<T,T_vector>
  {
  private:
    T_vector<T> cached_parents;
    T_vector<T> cached_value;
    std::mutex mtx;
  public:
    cached_dtm_node(size_t nparents,size_t ndim)
      :deterministic_node<T,T_vector>(nparents,ndim),cached_parents(nparents),cached_value(ndim)
    {}
    
    cached_dtm_node(size_t nparents)
      :deterministic_node<T,T_vector>(nparents,1),cached_parents(nparents),cached_value(1)
    {}
    
    cached_dtm_node()=delete;
    cached_dtm_node(const cached_dtm_node<T,T_vector>& )=delete;
    cached_dtm_node<T,T_vector>& operator=(const cached_dtm_node<T,T_vector>&)=delete;

    T do_value(size_t idx)const override
    {
      const_cast<std::mutex&>(mtx).lock();
      //T_vector<T> p(this->num_of_parents());
      T_vector<T> p(this->parent_values());
      bool parents_changed=false||(get_size(p)==0);
      
      for(size_t i=0;i<get_size(p);++i)
	{
	  ///set_element(p,i,this->parent(i));
	  if(get_element(p,i)!=get_element(cached_parents,i))
	    {
	      parents_changed=true;
	      
	      set_element(const_cast<cached_dtm_node<T,T_vector>*>(this)->cached_parents,i,get_element(p,i));
	      
	    }
	}
      
      
      T y=static_cast<T>(0);
      if(parents_changed)
	{
	  y=this->calc(idx,p);
	  set_element(const_cast<cached_dtm_node<T,T_vector>*>(this)->cached_value,idx,y);
	}
      else
	{
	  y=get_element(cached_value,idx);
	}
      const_cast<std::mutex&>(mtx).unlock();
      return y;
    }
  };
}


#endif
