#ifndef DUMP_GRAPH_TOPOLOGY
#define DUMP_GRAPH_TOPOLOGY
#include <core/graph.hpp>
#include <core/tag_t.hpp>
#include <iostream>
namespace mcmc_utilities
{
  template <typename T,template <typename TE> class T_vector>
  class topology_dumper
  {
  private:
    std::map<tag_t,T_vector<std::pair<tag_t,size_t> > > topology;
  public:
    topology_dumper(const graph<T,tag_t,T_vector>& g)
      :topology(g.topology())
    {}

    void all_to_dot(std::ostream& os)const
    {
      std::map<std::string,T_vector<std::string> > nodes;
      for(const auto& i:topology)
	{
	  std::string n1=i.first;
	  T_vector<std::string> plist;
	  for(const auto& j: i.second)
	    {
	      std::string n2=j.first;
	      push_back(plist,n2);
	    }
	  if(nodes.find(n1)==nodes.end())
	    {
	      nodes[n1]=plist;
	    }
	  
	}
      os<<"digraph G\n"
	<<"{\n";
      for(const auto& i:nodes)
	{
	  const T_vector<std::string> plist=i.second;
	  const std::string n=i.first;
	  if(plist.empty())
	    {
	      os<<n<<"\n";
	    }
	  else
	    {
	      for(const auto& j:plist)
		{
		  os<<j<<" -> "<<n<<"\n";
		}
	    }
	}
      os<<"}\n";
      
    }

    void to_dot(std::ostream& os)const
    {
      std::map<std::string,T_vector<std::string> > nodes;
      for(const auto& i:topology)
	{
	  std::string n1=i.first.name();
	  if(i.first.is_array())
	    {
	      n1+="_i";
	    }
	  T_vector<std::string> plist;
	  for(const auto& j: i.second)
	    {
	      std::string n2=j.first.name();
	      if(j.first.is_array())
		{
		  n2+="_i";
		}
	      push_back(plist,n2);
	    }
	  if(nodes.find(n1)==nodes.end())
	    {
	      nodes[n1]=plist;
	    }
	  
	}
      os<<"digraph G\n"
	<<"{\n";
      for(const auto& i:nodes)
	{
	  const T_vector<std::string> plist=i.second;
	  const std::string n=i.first;
	  if(plist.empty())
	    {
	      os<<n<<"\n";
	    }
	  else
	    {
	      for(const auto& j:plist)
		{
		  os<<j<<" -> "<<n<<"\n";
		}
	    }
	}
      os<<"}\n";
      
    }
  };

};

#endif
