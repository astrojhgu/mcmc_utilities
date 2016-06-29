#ifndef GRAPH_BUILDER
#define GRAPH_BUILDER
#include <memory>
#include <string>
#include <map>
#include <utility>
#include <cassert>
#include <set>

namespace mcmc_utilities
{
  
  
  template <typename T,template <typename T> class T_vector>
  class graph_builder
  {
  public:
    std::map<std::string,std::shared_ptr<vnode<T,T_vector> > > vnode_map;
    
  public:
    std::string add_node(const vnode<T,T_vector>& vn)
    {
      return add_node(vn.name,vn);
    }
    
    std::string add_node(std::string name,const vnode<T,T_vector>& vn)
    {
      auto p=vnode_map.find(name);
      if(p!=vnode_map.end())
	{
	  if(p->second->binded==true && vn.binded==true )
	    {
	      assert(0);
	    }
	}
      else if(vn.binded==true)
	{
	  auto p=vn.clone();
	  if(p->name!=name)
	    {
	      p->set_name(name);
	    }
	  vnode_map.insert(std::make_pair(name,p));
	}
      for(auto& i: vn.parents)
	{
	  if(i.first->binded)
	    {
	      add_node(i.first->name,*(i.first));
	    }
	}
      return name;
    }
    
    bool validate()
    {
      for(auto& p:vnode_map)
	{
	  for(auto& q:(p.second)->parents)
	    {
	      auto r=vnode_map.find(q.first->name);
	      if(r==vnode_map.end())
		{
		  //should never happen.
		  throw mcmc_exception("node not defined");
		}
	      assert((r->second->binded));
	      q.first=r->second;
	    }
	}
      return true;
    }

    void build(graph<T,std::string>& g)
    {
      for (auto& i: vnode_map)
	{
	  i.second->added=false;
	}
      for (auto& i: vnode_map)
	{
	  if (i.second->added==false)
	    {
	      add(g,i.second);
	    }
	}
    }

    void add(graph<T,std::string>& g,std::shared_ptr<vnode<T,T_vector> >& pn)
    {
      T_vector<std::pair<std::string,size_t> > parents;
      for(auto& i:pn->parents)
	{
	  if ( !(i.first->added))
	    {
	      add(g,i.first);
	    }
	  parents.push_back({i.first->name,i.second});
	}
      //cout<<pn->name<<" "<<typeid(*(pn.get())).name()<<" "<< pn->binded <<endl;
      g.add_node(pn->get_node(),pn->name,parents);
      pn->added=true;
    }
    
  };
  
  
  template <typename T>
  std::set<std::shared_ptr<vnode<T,T_vector> > > enumerate_all_named_parents(const graph_builder<T,T_vector>& gb,const std::shared_ptr<vnode<T,T_vector> >& pn)
  {
    std::set<std::shared_ptr<vnode<T,T_vector> > > result;
    for(auto & i : pn->parents)
      {
	if(i.first->named==true)
	  {
	    result.insert(i.first);
	  }
	else
	  {
	    for (auto& j:enumerate_all_named_parents(gb,i.first))

	      {
		result.insert(j);
	      }
	  }
      }
    return result;
  }

  std::string draw_node(std::string name,std::string type,int ninput,int noutput)
  {
    std::string result;
    result+=name;
    result+="[label=\"{";
    result+="{";
    for(int i=0;i<ninput;++i)
      {
	result+="<";
	result+="i";
	result+=std::to_string(i);
	result+=">";
	result+="i";
	result+=std::to_string(i);
	if(i!=ninput-1)
	  {
	    result+="|";
	  }
      }
    result+="}|";
    result+=name;
    result+=":";
    result+=type;
      
    result+="|{";
    for(int i=0;i<noutput;++i)
      {
	result+="<";
	result+="o";
	result+=std::to_string(i);
	result+=">";
	result+="o";
	result+=std::to_string(i);
	if(i!=noutput-1)
	  {
	    result+="|";
	  }
      }
    result+="}}\"];";
    return result;
  }
  
  template <typename T>
  void graph2dot1(const graph_builder<T,T_vector>& gb,std::ostream& os)
  {
    os<<"digraph{"<<std::endl;
    os<<"node [shape=record];"<<std::endl;
    for( const auto& i : gb.vnode_map )
      {
	std::string node_name = i.first;
	auto p=i.second->get_node();
	int nparents=p->num_of_parents();
	int ndim=p->num_of_dims();
	os<<draw_node(node_name,i.second->type,nparents,ndim)<<std::endl;;
      }
    
    for (const auto& i : gb.vnode_map)
      {
	std::string node1=i.first;
	int n=0;
	for (const auto& j: i.second->parents)
	  {
	    std::string node2=j.first->name;
	    os<<node2<<":o"<<j.second<<" -> "<<node1<<":i"<<n<<" ;"<<std::endl;
	    ++n;
	  }
      }
    
    os<<"}"<<std::endl;
  }
  
  template <typename T>
  void graph2dot2(const graph_builder<T,T_vector>& gb,std::ostream& os)
  {
    os<<"digraph{"<<std::endl;
    
    for (const auto& i : gb.vnode_map)
      {
	std::string node1=std::string("\"")+i.first+":"+i.second->type+"\"";
	if(i.second->named)
	  {
	    for (const auto& j:enumerate_all_named_parents(gb,i.second))
	      {
		std::string node2=std::string("\"")+j->name+":"+j->type+"\"";
		os<<node2<<" -> "<<node1<<" ;"<<std::endl;
		
	      }
	  }
	
      }
    
    os<<"}"<<std::endl;
  }
}


#endif

