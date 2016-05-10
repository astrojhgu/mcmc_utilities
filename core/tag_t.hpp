#ifndef TAG_T
#define TAG_T
#include <string>
#include <iostream>
namespace mcmc_utilities
{
  class tag_t
  {
  public:
    std::string _name;
    int _idx;
    bool _is_array;
  public:
    tag_t()
      :_name(),_idx(-1),_is_array(false)
    {}
    
    tag_t(std::string n)
      :_name(n),_idx(-1),_is_array(false)
    {}

    tag_t(std::string n,
	  int i)
      :_name(n),_idx(i),_is_array(true){}

    tag_t(std::string n,
	  int i,
	  bool ia)
      :_name(n),_idx(i),_is_array(ia){}

    tag_t(const tag_t& rhs)
      :_name(rhs._name),_idx(rhs._idx),_is_array(rhs._is_array)
    {
    }

    std::string name()const
    {
      return _name;
    }

    int idx()const
    {
      return _idx;
    }

    int is_array()const
    {
      return _is_array;
    }

    void set_name(const std::string& n)
    {
      _name=n;
    }

    void set_idx(int i)
    {
      _idx=i;
    }

    void set_is_array(bool b)
    {
      _is_array=b;
    }
    
    tag_t& operator=(const tag_t& rhs)
    {
      _name=rhs._name;
      _idx=rhs._idx;
      _is_array=rhs._is_array;
      return *this;
    }

    bool operator < (const tag_t& rhs)const
    {
      if(_name==rhs._name)
	{
	  return _idx<rhs._idx;
	}
      else
	{
	  return _name<rhs._name;
	}
    }

    bool operator == (const tag_t& rhs)const
    {
      if(!_is_array&&!rhs._is_array)
	{
	  return _name==rhs._name;
	}
      else
	{
	  return _name==rhs._name&&_idx==rhs._idx&&_is_array==rhs._is_array;
	}
    }

    bool operator != (const tag_t& rhs)const
    {
      return !(*this==rhs);
    }

    operator std::string()const
    {
      std::string result=_name;
      if(_is_array)
	{
	  result+="_";
	  result+=std::to_string(_idx);
	}
      return result;
    }
  };

  

  std::ostream& operator<<(std::ostream& os,const tag_t& t)
  {
    os<<t.name();
    if(t.is_array())
      {
	os<<"["<<t.idx()<<"]";
      }
    return os;
  }

};

#endif
