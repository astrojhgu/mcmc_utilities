#ifndef GRAPH_BUILDER
#define GRAPH_BUILDER
#include <memory>
#include <string>
#include <map>
#include <utility>
#include <cassert>


namespace mcmc_utilities
{
  
  
  template <typename T_p,typename T_var1>
  class graph_builder
  {
  public:
    std::map<std::string,shared_ptr<vnode<T_p,T_var1> > > vnode_map;
    
  public:
    void add_node(std::string name,const vnode<T_p,T_var1>& vn)
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
	  p->name=name;
	  vnode_map.insert(std::make_pair(name,p));
	}
      for(auto& i: vn.parents)
	{
	  if(i.first->binded)
	    {
	      add_node(i.first->name,*(i.first));
	    }
	}
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

    void build(graph<T_p,T_var1,std::string>& g)
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

    void add(graph<T_p,T_var1,std::string>& g,std::shared_ptr<vnode<T_p,T_var1> >& pn)
    {
      std::vector<std::pair<std::string,size_t> > parents;
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
  
  
  template <typename T_p,typename T_var1>
  void graph2gv(const graph_builder<T_p,T_var1>& gb,ostream& os)
  {
    os<<"digraph{"<<endl;
    
    for (const auto& i : gb.vnode_map)
      {
	std::string node1=string("\"")+i.first+":"+i.second->type+"\"";
	for (const auto& j: i.second->parents)
	  {
	    std::string node2=string("\"")+j.first->name+":"+j.first->type+"\"";
	    os<<node2<<" -> "<<node1<<" ;"<<endl;
	  }
      }
    
    os<<"}"<<endl;
  }
  
}


#endif

