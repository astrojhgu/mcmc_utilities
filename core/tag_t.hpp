#ifndef TAG_T
#define TAG_T
#include <string>
#include <iostream>
namespace mcmc_utilities
{
  struct tag_t
  {
    std::string name;
    int idx;
    bool is_array;
    
    tag_t()
      :name(),idx(-1),is_array(false)
    {}
    
    tag_t(std::string n)
      :name(n),idx(-1),is_array(false)
    {}

    tag_t(std::string n,
	  int i)
      :name(n),idx(i),is_array(true){}

    tag_t(std::string n,
	  int i,
	  bool ia)
      :name(n),idx(i),is_array(ia){}

    tag_t(const tag_t& rhs)
      :name(rhs.name),idx(rhs.idx),is_array(rhs.is_array)
    {
    }

    tag_t& operator=(const tag_t& rhs)
    {
      name=rhs.name;
      idx=rhs.idx;
      is_array=rhs.is_array;
      return *this;
    }

    bool operator < (const tag_t& rhs)const
    {
      if(name==rhs.name)
	{
	  return idx<rhs.idx;
	}
      else
	{
	  return name<rhs.name;
	}
    }

    bool operator == (const tag_t& rhs)const
    {
      if(!is_array&&!rhs.is_array)
	{
	  return name==rhs.name;
	}
      else
	{
	  return name==rhs.name&&idx==rhs.idx&&is_array==rhs.is_array;
	}
    }

    bool operator != (const tag_t& rhs)const
    {
      return !(*this==rhs);
    }

    operator std::string()const
    {
      std::string result=name;
      if(is_array)
	{
	  result+="_";
	  result+=std::to_string(idx);
	}
      return result;
    }
  };

  

  std::ostream& operator<<(std::ostream& os,const tag_t& t)
  {
    os<<t.name;
    if(t.is_array)
      {
	os<<"["<<t.idx<<"]";
      }
    return os;
  }

};

#endif
