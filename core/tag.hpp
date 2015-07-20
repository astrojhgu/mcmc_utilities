#ifndef TAG_HPP
#define TAG_HPP
#include <vector>

namespace mcmc_utilities
{
  template <typename T_str>
  class tag_t
  {
  private:
    T_str pname;
    std::vector<int> idx;
  public:
    tag_t()=delete;
    tag_t(const T_str& n,const std::vector<int>& idx_)
      :pname(n),idx(idx_)
    {}
    
    tag_t(const T_str& n)
      :pname(n),idx()
    {
      
    }

  public:
    T_str get_pname()const
    {
      return pname;
    }

    std::vector<int> get_idx()const
    {
      return idx;
    }
    
  public:
    bool operator<(const tag_t<T_str>& rhs)const
    {
      //return pname>=rhs.pname?true:
      if(pname==rhs.pname)
	{
	  if(idx.size()!=rhs.idx.size())
	    {
	      return idx.size()<rhs.idx.size();
	    }
	  for(int i=0;i<idx.size();++i)
	    {
	      if(idx[i]==rhs.idx[i])
		{
		  continue;
		}
	      return idx[i]<rhs.idx[i];
	    }
	}
      return (pname<rhs.pname);
    }
  };


  template <typename T_str>
  ostream& operator<<(ostream& os,const tag_t<T_str>& tag)
  {
    os<<tag.get_pname();
    std::vector<int> idx=tag.get_idx();
    for(int i=0;i<idx.size();++i)
      {
	os<<"("<<idx[i]<<")";
      }
    return os;
  }
}


#endif
