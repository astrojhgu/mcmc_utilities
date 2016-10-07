#ifndef WRAPPERED_GRAPH
#define WRAPPERED_GRAPH
#include <vector>

#include <core/graph.hpp>
#include <core/urand.hpp>

template <typename T>
using std_vector=std::vector<T>;

namespace mcmc_utilities
{
  template <typename T>
  class monitor_type
  {
  private:
    std::function<T()> func;
    
  public:
    monitor_type(const std::function<T()>& f)
    :func(f)
    {}
    
    T get()const
    {
      return func();
    }
  };
  
  template <typename T,typename T_tag>
  class wrappered_graph;
  
  template <typename T,typename T_tag>
  class node_adder
  {
  private:
    wrappered_graph<T,T_tag>& target;
    std_vector<std::pair<T_tag,size_t> > parents;
    std::shared_ptr<mcmc_utilities::node<double,std_vector> > pn;
    T_tag tag;
  public:
    node_adder(wrappered_graph<T,T_tag>& g,std::shared_ptr<mcmc_utilities::node<double,std_vector> > p,const T_tag& t)
      :target(g),pn(p),tag(t)
    {}

    
    node_adder& with_parent(const std::pair<T_tag,size_t>& p)
    {
      parents.push_back(p);
      return *this;
    }

    node_adder& with_parent(const T_tag& t,size_t n)
    {
      return with_parent(std::pair<T_tag,size_t>(t,n));
    }

    node_adder& with_parent(const T_tag& t)
    {
      return with_parent(std::pair<T_tag,size_t>(t,0));
    }
    

    node_adder& with_parent(const std::shared_ptr<mcmc_utilities::node<double,std_vector> >& p,size_t n)
    {
      T_tag t(target.get_tag(p));
      return with_parent(t,n);	      
    }

    node_adder& with_parent(const std::shared_ptr<mcmc_utilities::node<double,std_vector> >& p)
    {
      T_tag t(target.get_tag(p));
      return with_parent(t,0);
    }


    node_adder& with_parent(const std::pair<node<T,std_vector>*,size_t>& p)
    {
      return with_parent(std::shared_ptr<mcmc_utilities::node<double,std_vector> >(p.first),p.second);
    }

    node_adder& with_parent(node<T,std_vector>* p,size_t n)
    {
      return with_parent(std::shared_ptr<mcmc_utilities::node<double,std_vector> >(p),n);
    }

    node_adder& with_parent(node<T,std_vector>* p)
    {
      return with_parent(std::shared_ptr<mcmc_utilities::node<double,std_vector> >(p),0);
    }

    node_adder& with_value(const T& v)
    {
      auto p=std::dynamic_pointer_cast<stochastic_node<T,std_vector> >(pn);
      if(p)
	{
	  p->set_value(0,v);
	}
      return *this;
    }

    node_adder& with_value(const std_vector<T>& v)
    {
      auto p=std::dynamic_pointer_cast<stochastic_node<T,std_vector> >(pn);
      if(p)
	{
	  for(int i=0;i<v.size();++i)
	    {
	      p->set_value(i,v[i]);
	    }
	}
      return *this;
    }

    node_adder& with_observed_value(const T& v)
    {
      auto p=std::dynamic_pointer_cast<stochastic_node<T,std_vector> >(pn);
      if(p)
	{
	  p->set_observed(0,true);
	}
      return with_value(v);
    }

    node_adder& with_observed_value(const std_vector<T>& v)
    {
      auto p=std::dynamic_pointer_cast<stochastic_node<T,std_vector> >(pn);
      if(p)
	{
	  for(int i=0;i<v.size();++i)
	    {
	      p->set_value(i,v[i]);
	      p->set_observed(i,true);
	    }
	}
      return *this;
    }

    node_adder& with_value(size_t n,const T& v)
    {
      auto p=std::dynamic_pointer_cast<stochastic_node<T,std_vector> >(pn);
      if(p)
	{
	  p->set_value(n,v);
	}
      return *this;
    }

    node_adder& with_observed_value(size_t n,const T& v)
    {
      auto p=std::dynamic_pointer_cast<stochastic_node<T,std_vector> >(pn);
      if(p)
	{
	  p->set_observed(n,true);
	}
      return with_value(n,v);
    }

    void done()
      throw(node_name_already_used,parents_not_exist,invalid_node_type,
	    output_num_mismatch,parent_num_mismatch)
    {
      dynamic_cast<graph<T,T_tag,std_vector>&>(target).add_node(pn,tag,parents);
    }

    void with_tag(const T_tag& t)
      throw(node_name_already_used,parents_not_exist,invalid_node_type,
	    output_num_mismatch,parent_num_mismatch)
    {
      dynamic_cast<graph<T,T_tag,std_vector>&>(target).add_node(pn,t,parents);
    }
  };
  
  template <typename T,typename T_tag>
  class wrappered_graph
    :public graph<T,T_tag,std_vector>
  {
  public:
    urand<T> rnd;
  public:
    void sample()
    {
      graph<T,T_tag,std_vector>::sample(rnd);
    }

    monitor_type<T> get_monitor(const T_tag& t,size_t n)
    {
      return monitor_type<T>(graph<T,T_tag,std_vector>::get_monitor(t,n));
    }

    T get_value(const T_tag& t,size_t n)
    {
      return graph<T,T_tag,std_vector>::get_value(t,n);
    }

    void set_value(const T_tag& t,size_t n,T v)
    {
      graph<T,T_tag,std_vector>::set_value(t,n,v);
    }

    bool is_observed(const T_tag& t,size_t n)
    {
      return graph<T,T_tag,std_vector>::is_observed(t,n);
    }

    void set_observed(const T_tag& t,size_t n,bool b)
    {
      graph<T,T_tag,std_vector>::set_observed(t,n,b);
    }

    node_adder<T,T_tag> add_node(std::shared_ptr<mcmc_utilities::node<double,std_vector> >& p,const T_tag& t=T_tag())
      throw(node_name_already_used,parents_not_exist,invalid_node_type,
	    output_num_mismatch,parent_num_mismatch)
    {
      return node_adder<T,T_tag>(*this,p,t);
    }

    node_adder<T,T_tag> add_node(mcmc_utilities::node<T,std_vector>* p,const T_tag& t=T_tag())
      throw(node_name_already_used,parents_not_exist,invalid_node_type,
	    output_num_mismatch,parent_num_mismatch)
    {
      return node_adder<T,T_tag>(*this,std::shared_ptr<mcmc_utilities::node<double,std_vector> >(p),t);
    }
    
  };

  
}

#endif
