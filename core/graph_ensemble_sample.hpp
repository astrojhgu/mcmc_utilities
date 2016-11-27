#ifndef GRAPH_ENSEMBLE_SAMPLE
#define GRAPH_ENSEMBLE_SAMPLE
#include <memory>
#include "ensemble_sample.hpp"
namespace mcmc_utilities
{
  template <typename T,typename T_tag,template <typename TE> class T_vector>
  class graph_function_adapter
  {
  private:
    std::shared_ptr<graph<T,T_tag,T_vector> > g;
    std::vector<std::shared_ptr<graph_function_adapter<T,T_tag,T_vector> > > copies;
    size_t fetch_count;
  public:
    graph_function_adapter(const graph<T,T_tag,T_vector>& g1,size_t nthread)
      :g(new graph<T,T_tag,T_vector>()),copies(),fetch_count(0)
    {
      g->copy_from(g1);
      if(nthread>1)
	{
	  for(size_t i=0;i<nthread;++i)
	    {
	      auto g2=new graph_function_adapter<T,T_tag,T_vector>(g1,1);
	      copies.push_back(std::shared_ptr<graph_function_adapter<T,T_tag,T_vector> >(g2));
	    }
	}
    }

    graph_function_adapter(const graph_function_adapter<T,T_tag,T_vector>&)=default;
    graph_function_adapter<T,T_tag,T_vector>& operator=(const graph_function_adapter<T,T_tag,T_vector>&)=default;

    graph_function_adapter<T,T_tag,T_vector>& fetch_copy()
    {
      if(copies.empty())
	{
	  mcmc_exception e("should never called when copies is empty");
	  throw e;
	}
      fetch_count++;
      fetch_count%=copies.size();
      return *(copies[fetch_count]);
    }

    T operator()(const T_vector<T>& x)
    {
      return g->eval_logprob(x);
    }
  };

  template <typename T,typename T_tag,template <typename TE> class T_vector>
  graph_function_adapter<T,T_tag,T_vector>& clone(graph_function_adapter<T,T_tag,T_vector>& gfa)
  {
    return gfa.fetch_copy();
  }


  
  template <typename T,typename T_graph,typename T_ensemble>
  T_ensemble ensemble_sample(T_graph& g,
			      const T_ensemble& ensemble,
			      base_urand<T>& rng)
  {
    typedef typename element_type_trait<T_ensemble>::element_type T_var;
    auto p=g.get_params();
    size_t nparams=get_size(p);
    auto ensemble1=clone(ensemble);
    while(get_size(ensemble1)<2*nparams)
      {
	g.sample(rng);
	auto p=g.get_params();
	if(get_size(ensemble1)>0)
	  {
	    auto plast=get_element(ensemble1,get_size(ensemble1)-1);
	    bool no_jump=true;
	    for(size_t i=0;i<nparams;++i)
	      {
		if(get_element(p,i)!=get_element(plast,i))
		  {
		    no_jump=false;
		    break;
		  }
	      }
	    if(no_jump)
	      {
		continue;
	      }
	  }
	push_back(ensemble1,p);
      }
    ensemble1=ensemble_sample([&g](const T_var& x)
			     {
			       T result=g.eval_logprob(x);
			       return result;
			     },ensemble1,rng);
    return ensemble1;
  }
}


#endif
