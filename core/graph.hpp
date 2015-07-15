#include <memory>
#include <vector>
#include <map>
#include "mcmc_exception.hpp"
#include "gibbs_sampler.hpp"

namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class stochastic_node;

  template <typename T_p,typename T_var1>
  class deterministic_node;
  
  template <typename T_p,typename T_var1>
  class node
  {
    friend class stochastic_node<T_p,T_var1>;
    friend class deterministic_node<T_p,T_var1>;
  protected:
    std::vector<stochastic_node<T_p,T_var1>* > stochastic_children;
    std::vector<deterministic_node<T_p,T_var1 >* > deterministic_children;
    std::vector<node<T_p,T_var1>* > parents;

  public:
    explicit node(int nparents)
      :parents(nparents)
    {}

    node()=delete;
    node(const node<T_p,T_var1>& rhs)=delete;
    node<T_p,T_var1>& operator=(const node<T_p,T_var1>& rhs)=delete;

    virtual ~node(){}
    
  public:
    int num_of_parents()const
    {
      return parents.size();
    }
    
    T_var1 value()const
    {
      return do_value();
    }

    virtual T_p log_likelihood()const final
    {
      T_p result=0;
      for(auto& p : stochastic_children)
	{
	  result+=p->log_prior_prob();
	}
      for(auto& p:deterministic_children)
	{
	  result+=p->log_likelihood();
	}
      return result;
    }

    void connect_to_parent(node<T_p,T_var1>* prhs,int n)
    {
      do_connect_to_parent(prhs,n);
    }

    
    
  private:
    virtual T_var1 do_value()const=0;
    virtual void do_connect_to_parent(node<T_p,T_var1>* prhs,int n)=0;
  };

  
  template <typename T_p>
  class base_urand
  {
  public:
    T_p operator()()const
    {
      return do_rand();
    }
  private:
    virtual T_p do_rand()const=0;
      
  };
  

  template <typename T_p,typename T_var1>
  class stochastic_node
    :public node<T_p,T_var1>,public probability_density_1d<T_p,T_var1>
  {
  private:
    T_p v;

  public:
    stochastic_node(int nparents,T_var1 v_)
      :node<T_p,T_var1>(nparents),
      v(v_)
    {}

    stochastic_node()=delete;
    stochastic_node(const stochastic_node<T_p,T_var1>& )=delete;
    stochastic_node<T_p,T_var1>& operator=(const stochastic_node<T_p,T_var1>&)=delete;

  public:
    T_p log_post_prob()const
    {
      return log_prior_prob()+this->log_likelihood();
    }

    T_p do_eval_log(const T_var1& x)const
    {
      const_cast<stochastic_node*>(this)->set_value(x);
      return log_post_prob();
    }
    
    T_p log_prior_prob()const
    {
      return do_log_prior_prob();
    }

    void set_value(const T_var1& v_)
    {
      v=v_;
    }

  public:
    void sample(const base_urand<T_p>& rnd)
    {
      do_sample(rnd);
    }
  private:
    virtual T_p do_log_prior_prob()const=0;
    T_var1 do_value()const override
    {
      return v;
    }

    void do_connect_to_parent(node<T_p,T_var1>*  rhs,int n) override
    {
      this->parents.at(n)=rhs;
      rhs->stochastic_children.push_back(this);
    }

    virtual void do_sample(const base_urand<T_p>&)=0;
  };

  template <typename T_p,typename T_var1>
  class continuous_node
    :public stochastic_node<T_p,T_var1>
  {
  public:
    continuous_node(int nparents,T_var1 v_)
      :stochastic_node<T_p,T_var1>(nparents,v_)
    {}
    
    continuous_node()=delete;
    continuous_node(const continuous_node<T_p,T_var1>& )=delete;
    continuous_node<T_p,T_var1>& operator=(const continuous_node<T_p,T_var1>&)=delete;

  public:
    void do_sample(const base_urand<T_p>& rnd)override
    {
      T_var1 xprev=this->value();
      constexpr int nsamp=10;
      std::vector<T_var1> xsamp(nsamp);
      arms_simple(*this,xprev,xsamp,dometrop(),rnd);
      this->set_value(xsamp.back());
    }
  private:
    virtual bool dometrop()const
    {
      return true;
    }
  };


  template <typename T_p,typename T_var1>
  class discrete_node
    :public stochastic_node<T_p,T_var1>
  {
  public:
    discrete_node(int nparents,T_var1 v_)
      :stochastic_node<T_p,T_var1>(nparents,v_)
    {}
    
    discrete_node()=delete;
    discrete_node(const discrete_node<T_p,T_var1>& )=delete;
    discrete_node<T_p,T_var1>& operator=(const discrete_node<T_p,T_var1>&)=delete;

  public:
    void do_sample(const base_urand<T_p>& rnd)override
    {
      this->set_value(discrete_sample(*this,rnd));
    }
  };

  
  
  template <typename T_p,typename T_var1>
  class deterministic_node
    :public node<T_p,T_var1>
  {
  public:
    deterministic_node(int nparents)
      :node<T_p,T_var1>(nparents)
    {}
    
    deterministic_node()=delete;
    deterministic_node(const deterministic_node<T_p,T_var1>& )=delete;
    deterministic_node<T_p,T_var1>& operator=(const deterministic_node<T_p,T_var1>&)=delete;

  private:
    void do_connect_to_parent(node<T_p,T_var1>*  rhs,int n) override
    {
      this->parents.at(n)=rhs;
      rhs->deterministic_children.push_back(this);
    }
  };

  
  template <typename T_p,typename T_var1,typename T_str>
  class graph
  {
  private:
    std::vector<stochastic_node<T_p,T_var1>* > stochastic_node_list;
    std::vector<deterministic_node<T_p,T_var1>* > deterministic_node_list;
    std::vector<node<T_p,T_var1>* > observed_node_list;
    std::map<std::pair<T_str,int>,std::shared_ptr<node<T_p,T_var1> > > node_map;

  public:
    void clear()
    {
      node_map.clear();
      stochastic_node_list.clear();
      deterministic_node_list.clear();
      observed_node_list.clear();
    }

    void sample(const base_urand<T_p>& urand)
    {
      for(auto& p:stochastic_node_list)
	{
	  p->sample(urand);
	}
    }

    std::vector<T_var1> get_params()const
    {
      std::vector<T_var1> result(stochastic_node_list.size());
      for(int i=0;i<result.size();++i)
	{
	  result[i]=stochastic_node_list[i]->value();
	}
      return result;
    }

    void set_value(const std::pair<T_str,int>& tag,
		      const T_var1& v)
    {
      if(node_map.count(tag)==0)
	{
	  throw mcmc_exception("node not found by tag");
	}
      
      stochastic_node<T_p,T_var1>* ps=dynamic_cast<stochastic_node<T_p,T_var1>* >(node_map[tag].get());
      if(ps==nullptr)
	{
	  throw mcmc_exception("this node is not a stochastic node");
	}
      ps->set_value(v);

    }

    T_p log_likelihood(const std::pair<T_str,int>& tag)const
    {
      auto i=node_map.find(tag);
      if(i==node_map.end())
	{
	  throw mcmc_exception("node not found");
	}
      i->second->log_likelihood();
    }

    T_p log_post_prob(const std::pair<T_str,int>& tag)const
    {
      auto i=node_map.find(tag);
      if(i==node_map.end())
	{
	  throw mcmc_exception("node not found");
	}

      auto ps=dynamic_cast<stochastic_node<T_p,T_var1>*>(i->second.get());
      if(ps==nullptr)
	{
	  throw mcmc_exception("this node is not a stochastic node");
	}
      return ps->log_post_prob();
    }

        
    void add_node(node<T_p,T_var1>* pn,bool observed,
		  const std::pair<T_str,int>& tag,
		  const std::vector<std::pair<T_str,int> >& parents)
    {
      std::shared_ptr<node<T_p,T_var1> > ptr(pn);
      if (node_map.count(tag)!=0)
	{
	  throw mcmc_exception("node name already exists");
	}
      if(pn->num_of_parents()!=parents.size())
	{
	  throw mcmc_exception("parent num do not match");
	}
      for(auto& p:parents)
	{
	  if(node_map.count(p)==0)
	    {
	      throw mcmc_exception("parent should be added before hand");
	    }
	}
      
      auto ps=dynamic_cast<stochastic_node<T_p,T_var1>*>(pn);
      if(ps!=nullptr)
	{
	  if(!observed)
	    {
	      stochastic_node_list.push_back(ps);
	    }
	  else
	    {
	      observed_node_list.push_back(pn);
	    }
	}
      else
	{
	  auto pd=dynamic_cast<deterministic_node<T_p,T_var1>* >(pn);
	  if(pd==nullptr)
	    {
	      throw mcmc_exception("input node is neither stochastic node, nor deterministic node");
	    }
	  deterministic_node_list.push_back(pd);
	}
      int n=0;
      for(auto& i:parents)
	{
	  auto n_iter=node_map.find(i);
	  pn->connect_to_parent(n_iter->second.get(),n);
	  ++n;
	}
      node_map[tag]=ptr;
    }
  };
}
