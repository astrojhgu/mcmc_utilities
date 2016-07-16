#ifndef TP_AWARE_NODE_HPP
#define TP_AWARE_NODE_HPP
#include "cached_dtm_node.hpp"
#include <set>
#include <stack>

namespace mcmc_utilities
{
  struct order
  {
    int n{0};
    bool homo{true};
    bool poly{true};
    order(int _n,bool _homo,bool _poly)
      :n(_n),
       homo(_homo),
       poly(_poly)
    {}

  public:
    int get_n()const
    {
      return n;
    }

    bool is_homo()const
    {
      return homo;
    }

    bool is_poly()const
    {
      return poly;
    }
  };

  std::ostream& operator<<(std::ostream& os,const order& o)
  {
    os<<"("
      <<o.n
      <<","
      <<(o.is_homo()?std::string("homo"):std::string("non-homo"))
      <<","
      <<(o.is_poly()?std::string("poly"):std::string("non-poly"))<<")";
    return os;
  }

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


  public:
    order get_order(const node<T,T_vector>* pn,int n)const
    {
      return do_get_order(pn,n);
    }

    order get_parent_order(int m,const node<T,T_vector>* pn,int n)const
    {
      //auto p=get_element(this->parents,m).first;
      auto p=this->get_parent(m).first;
      if(p==pn)
	{
	  return order{this->get_parent(m).second==n?1:0,true,true};
	}
      
      auto x= dynamic_cast<const tp_aware_dtm_node<T,T_vector>*>(this->get_parent(m).first);
      if(x==nullptr)
	{
	  return order{0,true,true};
	}
      else
	{
	  return x->get_order(pn,n);
	}
    }
    
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

  private:
    virtual order do_get_order(const node<T,T_vector>* pn,int n)const=0;
  };
}


#endif
