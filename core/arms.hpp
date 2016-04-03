#ifndef MYARMS_HPP
#define MYARMS_HPP
//#include "distribution.hpp"
#include "error_handler.hpp"
#include <algorithm>
#include <limits>
#include <iomanip>
#include <iostream>
//#define DEBUG
#ifdef DEBUG
#include <fstream>
#include <logger/logger.hpp>
#endif


#include <cassert>
#include <list>
#include <utility>
#include <limits>
#include <mutex>
//all the y in following program actually equal to log p(x)

namespace mcmc_utilities
{
  template <typename T,typename TD>
  T eval_log(const TD& pd,T x,T scale)
  {
    T result=pd(x)-scale;
#ifdef DEBUG
    if(std::isinf(result)||std::isnan(result))
      {
	std::cerr<<result<<std::endl;
	std::cerr<<scale<<std::endl;
	std::cerr<<(&pd)<<std::endl;
	std::cerr<<pd(x)<<std::endl;
	assert(0);
      }
    
#endif
    if(std::isnan(result)||std::isinf(result))
      {
	nan_or_inf e;
	e.attach_message("#0");
	throw e;
      }
    return result;
  }
  
  template <typename T>
  T int_exp_y(const T& x,const std::pair<T,T>& p1,const std::pair<T,T>& p2)
  {
    
    //calculate \int_x1^2 exp(y1+(x-x1)*(y2-y1)/(x2-x1)) dx
    //i.e., cdf of distribution e^h(x)
    const T& x1=p1.first;
    const T& x2=p2.first;
    const T& y1=p1.second;
    const T& y2=p2.second;
    //std::cerr<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<std::endl;
    T result=0;
    if(x1==x2)
      {
	//result=0;
      }
    else if(!std::isinf(y1)&&!std::isinf(y2))
      {
	const T k=(y2-y1)/(x2-x1);
	if(k==0)
	  {
	    result=std::exp(y1)*(x-x1);
	  }
	else
	  {
	    result=(std::exp(k*(x-x1))-1)*std::exp(y1)/k;
	    if(std::isnan(result))
	      {
		result=(std::exp(k*(x-x1)+y1)-std::exp(y1))/k;
	      }
	    
	  }
#ifdef DEBUG
	if(std::isnan(result))
	  {
	    std::cerr<<std::setprecision(20)<<k<<" "<<x1<<" "<<y1<<" "<<x<<std::endl;
	    std::cerr<<std::exp(k*(x-x1)+y1)<<std::endl;
	    std::cerr<<(std::exp(k*(x-x1)+y1)-std::exp(y1))<<std::endl;
	    std::cerr<<(std::exp(k*(x-x1))-1)*std::exp(y1)<<std::endl;
	    assert(0);
	  }
#endif
      }
    else
      {
	result=0;
      }
    
    if(std::isnan(result))
      {
	//nan_or_inf e;
	//e.attach_message("nan #1");
	ill_conditioned_distribution e;
	e.attach_message("ill #1");
	throw e;
      }
    return result;
    /*
    if(y1!=y2)
      {
	T result=0;
	const T a=(x*y1-x*y2+x1*y2-x2*y1)/(x1-x2)-y1;
	if(a<static_cast<T>(30))
	  {
	    result=(std::exp(a)-1)*(x1-x2)/(y1-y2)*std::exp(y1);
	  }
	else
	  {
	    result=std::exp(a+y1)*(x1-x2)/(y1-y2);
	  }
	if(std::isnan(result))
	  {
	    std::cerr<<"x1="<<x1<<" y1="<<y1<<" x2="<<x2<<" y2="<<y2<<" a="<<a<<std::endl;
	    assert(0);
	    throw nan_or_inf();
	  }
	result=std::max(static_cast<T>(0),result);
	assert(!std::isnan(result));
	//assert(!std::isinf(result));
	return result;
      }
    else
      {
	T result=std::exp(y1)*(x-x1);
	assert(!std::isnan(result));
	//assert(!std::isinf(result));
	return result;
      }
    */
  }

  template <typename T>
  T inv_int_exp_y(const T& Z,const std::pair<T,T>& p1,const std::pair<T,T>& p2)
  {
    //
    //solve[\int_x1^x exp(k*(x-x1)+y1) dx==Z,x]

    const T& x1=p1.first;
    const T& x2=p2.first;
    const T& y1=p1.second;
    const T& y2=p2.second;
    
    T result=0;
    T k=(y2-y1)/(x2-x1);
    if(x1==x2||Z==0)
      {
	result=x1;
      }
    else
      {	
	if(k==0)
	  {
	    result= x1+Z*std::exp(-y1);
	    if(std::isinf(result))
	      {
		result=x1+std::exp(std::log(Z)-y1);
	      }
	  }
	else
	  {
	    T U=1+k*Z*std::exp(-y1);
	    
	    
	    if(!std::isinf(U))
	      {
		T logU=0;
		if(U<=0)
		  {
		    logU=std::log((std::exp(y1)+k*Z))-y1;
		  }
		else
		  {
		    logU=std::log(U);
		  }
		result=x1+logU/k;
	      }
	    else
	      {
		result=x1+(std::log(k*Z)-y1)/k;
	      }
	  }
      }
#ifdef DEBUG
    if(std::isinf(result)||std::isnan(result))
      {
	std::cerr<<std::setprecision(20);
	std::cerr<<"x1="<<x1<<" y1="<<y1<<std::endl;
	std::cerr<<"x2="<<x2<<" y2="<<y2<<std::endl;
	std::cerr<<"Z="<<Z<<std::endl;
	std::cerr<<k*Z<<std::endl;
	std::cerr<<y1<<std::endl;
	std::cerr<<1+k*Z*std::exp(-y1)<<std::endl;
	assert(0);

	//throw nan_or_inf();
      }
#endif

    if(!std::isfinite(result))
      {
	nan_or_inf e;
	e.attach_message("#2");
	throw e;
      }
    
    return result;

    /*
    if(y2!=y1&&x1!=x2)
      {
	//result=x1 + (x1 - x2)/(y1 - y2) * std::log( (1 + (y1 - y2)/(x1 - x2)*y *std::exp(-y1) ));
	
	T result=0;
	T f1=0;
	T f2=0;
	
	//T f1=std::log(1 + (y2 - y1)/(x2 - x1)*y *std::exp(-y1));
	                ~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^f2
			
	f2=(x2-x1+y *std::exp(-y1) * (y2 - y1))/(x2 - x1);
	
	
	//T f3=(y2>y1)?(1+std::exp(f2)):(1-std::exp(f2));
	if(!isinf(f2))
	  {
	    f1=std::log(f2);
	  }
	else if(f2>0)
	  {
	    f1=std::log(y2-y1)-std::log(x2-x1)+std::log(y)-y1;
	  }
	else
	  {
	    assert(0);
	  }
	
	
	//result=x1 + (x1 - x2)/(y1 - y2) * std::log(1 + (y1 - y2)/(x1 - x2)*y *std::exp(-y1));
	result=x1 + (x1 - x2)/(y1 - y2) * f1;
	if(std::isnan(result))
	  {
	    std::cerr<<std::setprecision(20)<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<" "<<y<<std::endl;
	    assert(0);
	    throw nan_or_inf();
	  }
	if(std::isinf(result))
	  {
	    std::cerr<<std::setprecision(20)<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<" "<<y<<std::endl;
	    assert(0);
	    throw nan_or_inf();
	  }
	return result;
      }
    else if(y1==y2&&x1!=x2)
      {
	return y*std::exp(-y1)+x1;
      }
      return x1;
    */
  }
  
  template <typename T>
  class section
  {
  private:
    T _x_l;//lower x of this section
    T _x_u;//upper x of this section
    T _x_i;//intersection x of this section
    
    T _y_l;//ln y(x_l) of this section
    T _y_u;//ln y(x_u) of this section
    T _y_i;//ln y(x_i) of this section

    T _int_exp_y_l;//\int_x_l^x_i exp(y)
    T _int_exp_y_u;//\int_x_i^x_u exp(y)
    T _cum_int_exp_y_l;
    T _cum_int_exp_y_u;

  public:
    T x_l()const
    {
      return _x_l;
    }

    T x_u()const
    {
      return _x_u;
    }

    T x_i()const
    {
      return _x_i;
    }
   
    T y_l()const
    {
      return _y_l;
    }

    T y_u()const
    {
      return _y_u;
    }

    T y_i()const
    {
      return _y_i;
    } 
    
    T int_exp_y_l()const
    {
      return _int_exp_y_l;
    }

    T int_exp_y_u()const
    {
      return _int_exp_y_u;
    }

    T cum_int_exp_y_l()const
    {
      return _cum_int_exp_y_l;
    }
    
    T cum_int_exp_y_u()const
    {
      return _cum_int_exp_y_u;
    }

    //////////
    void set_x_l(T x)
    {
      _x_l=x;
    }

    void set_x_u(T x)
    {
      _x_u=x;
    }

    void set_x_i(T x)
    {
#ifdef DEBUG
      assert(!std::isnan(_x_l));
      assert(!std::isnan(_x_u));
#endif      
      _x_i=x;
      if(_x_i<_x_l)
	{
	  _x_i=_x_l;
	}
      if(_x_i>_x_u)
	{
	  _x_i=_x_u;
	}
    }
   
    void set_y_l(T y)
    {
#ifdef DEBUG
      assert(!std::isnan(y));
#endif
      _y_l=y;
    }

    void set_y_u(T y)
    {
#ifdef DEBUG
      assert(!std::isnan(y));
#endif
      _y_u=y;
    }

    void set_y_i(T y)
    {
#ifdef DEBUG
      assert(!std::isnan(y));
#endif
      _y_i=y;
    } 
    
    void set_int_exp_y_l(T y)
    {
#ifdef DEBUG
      assert(!std::isnan(y));
#endif
      _int_exp_y_l=y;
    }

    void set_int_exp_y_u(T y)
    {
#ifdef DEBUG
      assert(!std::isnan(y));
#endif

      _int_exp_y_u=y;
    }

    void set_cum_int_exp_y_l(T y)
    {
#ifdef DEBUG
      assert(!std::isnan(y));
#endif

      _cum_int_exp_y_l=y;
    }
    
    void set_cum_int_exp_y_u(T y)
    {
#ifdef DEBUG
      assert(!std::isnan(y));
#endif
      _cum_int_exp_y_u=y;
    }
    
  public:
    section()
      :_x_l(std::nan("")),_x_u(std::nan("")),_x_i(std::nan("")),
       _y_l(std::nan("")),_y_u(std::nan("")),_y_i(std::nan("")),
       _int_exp_y_l(std::nan("")),
       _int_exp_y_u(std::nan("")),
       _cum_int_exp_y_l(std::nan("")),
       _cum_int_exp_y_u(std::nan(""))
    {}
    
  public:
    bool encloses(const T& x)const
    {
      return x>=_x_l&&x<=_x_u;
    }
    
    void calc_int_exp_y()
    {
      if(_x_i==_x_l)
	{
	  _int_exp_y_l=0;
	}
      else
	{
#ifdef DEBUG
	  assert(!std::isinf(_y_i)||_y_i<0);
	  assert(!std::isnan(_y_i));
#endif
	  
	  _int_exp_y_l=int_exp_y(_x_i,std::make_pair(_x_l,_y_l),std::make_pair(_x_i,_y_i));
	  
#ifdef DEBUG
	  assert(!std::isnan(_int_exp_y_l));
#endif
	  //assert(!std::isinf(int_exp_y_l));
	  
	}
      if(_x_i==_x_u)
	{
	  _int_exp_y_u=0;
	}
      else
	{
#ifdef DEBUG
	  assert(!std::isinf(_y_i)||_y_i<0);
	  assert(!std::isnan(_y_i));
#endif
	  
	  _int_exp_y_u=int_exp_y(_x_u,std::make_pair(_x_i,_y_i),std::make_pair(_x_u,_y_u));
	  
#ifdef DEBUG
	  assert(!std::isnan(_int_exp_y_u));
#endif
	  //assert(!std::isinf(int_exp_y_u));
	  
	}
    }

  public:
    T eval_y(const T& x)const
    {
      if(x==_x_l)
	{
	  return _y_l;
	}
      if(x==_x_u)
	{
	  return _y_u;
	}
      if(x==_x_i)
	{
	  return _y_i;
	}
      if(x<_x_i)
	{
	  return _y_l+(x-_x_l)*(_y_i-_y_l)/(_x_i-_x_l);
	}
      else
	{
	  return _y_u+(x-_x_u)*(_y_i-_y_u)/(_x_i-_x_u);
	}
    }
  };

  template <typename T>
  std::ostream& operator<<(std::ostream& os,const section<T>& s)
  {
    //os<<s.x_i()<<" "<<s.cum_int_exp_y_l()<<" "<<s.y_i()<<"\n";
    //os<<s.x_u()<<" "<<s.cum_int_exp_y_u()<<" "<<s.y_u()<<"\n";
    os<<"!\n";
    os<<s.x_l()<<" "<<s.y_l()<<"\n";
    os<<s.x_i()<<" "<<s.y_i()<<"  !"<<s.int_exp_y_l()<<"\n";
    os<<s.x_u()<<" "<<s.y_u()<<"  !"<<s.int_exp_y_u()<<"\n";

    return os;
  }
  
  template <typename T>
  std::ostream& operator<<(std::ostream& os,const std::list<section<T> >& section_list)
  {
    if(section_list.empty())
      {
	return os;
      }
    os<<section_list.front().x_l()<<" "<<0<<" "<<section_list.front().y_l()<<"\n";
    for(auto& i:section_list)
      {
	os<<i;
      }
    return os;
  }
  

  template <typename T>
  T eval(const T& x,const std::list<section<T> >& section_list)
  {
    for(auto& i:section_list)
      {
	if(i.encloses(x))
	  {
	    return i.eval_y(x);
	  }
      }
    //cerr<<"x="<<x<<endl;
    var_out_of_range e;
    e.attach_message("#3");
    throw e;
    return 0;
  }

  template <typename T>
  T eval_ey(const T& x,const std::list<section<T> >& section_list)
  {
    return std::exp(eval(x,section_list));
  }
  
  template <typename T>
  std::pair<T,T> solve_intersection(const section<T>& s,
				    const std::pair<T,T>& p1,
				    const std::pair<T,T>& p2,
				    const std::pair<T,T>& p3,
				    const std::pair<T,T>& p4)
  {
    const T& x1=p1.first;
    const T& x2=p2.first;
    const T& x3=p3.first;
    const T& x4=p4.first;
    
    const T& y1=p1.second;
    const T& y2=p2.second;
    const T& y3=p3.second;
    const T& y4=p4.second;

#ifdef DEBUG
    assert(!std::isnan(y1)&&!std::isnan(y2)&&!std::isnan(y3)&&!std::isnan(y4));
#endif
    

    T x_i,y_i;

    if(std::isinf(y2))
      {
	x_i=x1;
	y_i=y4+(x1-x4)/(x3-x4)*(y3-y4);

	if(y_i<y1)
	  {
	    x_i=(s.x_l()+s.x_u())/2;
	    y_i=(s.y_l()+s.y_u())/2;
	  }
	//assert(!std::isinf(y_i));
      }
    else if(std::isinf(y4))
      {
	x_i=x3;
	y_i=y1+(x3-x1)/(x2-x1)*(y2-y1);

	if(y_i<y3)
	  {
	    x_i=(s.x_l()+s.x_u())/2;
	    y_i=(s.y_l()+s.y_u())/2;
	  }
	//assert(!std::isinf(y_i));
      }
    else
      {

    
	x_i=-((-(x3 - x4) * (x2 * y1 - x1 * y2) + (x1 - x2) * (x4 * y3 -  x3 * y4))/(-(x3 - x4) * (-y1 + y2) + (x1 - x2) * (-y3 + y4)));
	y_i=-(x2 * y1 * y3 - x4 * y1 * y3 - x1 * y2 * y3 + x4 * y2 * y3 - x2 * y1 * y4 + x3 * y1 * y4 +  x1 * y2 * y4 - x3 * y2 * y4)/(-x3 * y1 + x4 * y1 + x3 * y2 - x4 * y2 + x1 * y3 - x2 * y3 - x1 * y4 + x2 * y4);

	if(((y3-y1)*(x2-x1)==(y2-y1)*(x3-x1)&&
	    (y4-y1)*(x2-x1)==(y2-y1)*(x4-x1))||(y2-y1)*(x4-x3)==(x2-x1)*(y4-y3)||!std::isfinite(y_i))
	  {
	    x_i=(s.x_l()+s.x_u())/2;
	    y_i=(s.y_l()+s.y_u())/2;
	  }
	
	if((x_i<=s.x_l())||(x_i>=s.x_u())||y_i<(y2+(x_i-x2)/(x3-x2)*(y3-y2)))
	  {
	    x_i=(s.x_l()+s.x_u())/2;
	    y_i=(s.y_l()+s.y_u())/2;
	  }
#ifdef DEBUG
	assert(!std::isinf(y_i)||y_i<0);
#endif
      }
    
    if(std::isnan(x_i)||std::isnan(y_i))
      {
#ifdef DEBUG	
	std::cerr<<std::setprecision(20)<<x1<<" "<<y1<<std::endl;
	std::cerr<<x2<<" "<<y2<<std::endl;
	std::cerr<<x3<<" "<<y3<<std::endl;
	std::cerr<<x4<<" "<<y4<<std::endl;
	//assert(0);
#endif
	nan_or_inf e;
	e.attach_message("#4");
	throw e;
	    
      }
    return std::make_pair(x_i,y_i);
  }

  template <typename T>
  void calc_cum_int_exp_y(std::list<section<T> >& section_list,typename std::list<section<T> >::iterator& i)
  {
    i->calc_int_exp_y();
    T cum_from=0;
    if(i!=section_list.begin())
      {
	auto i_prev=i;
	i_prev--;
	cum_from=i_prev->cum_int_exp_y_u();
      }
    i->set_cum_int_exp_y_l(cum_from+i->int_exp_y_l());
    i->set_cum_int_exp_y_u(i->cum_int_exp_y_l()+i->int_exp_y_u());

#ifdef DEBUG
    assert(!std::isnan(i->cum_int_exp_y_u()));
    assert(i->cum_int_exp_y_u()>=0);
#endif

    if(std::isnan(i->cum_int_exp_y_u()))
      {
#ifdef DEBUG
	assert(0);
#endif
	nan_or_inf e;
	e.attach_message("#5");
	throw e;
      }

    if(i->cum_int_exp_y_u()<0)
      {
	cum_lt_zero e;
	e.attach_message("#6");
	throw e;
      }
  }

  template <typename T>
  void calc_intersection(std::list<section<T> >& section_list,typename std::list<section<T> >::iterator& i)
  {
    auto i_next=i;
    auto i_prev=i;
    
    if(i_next!=section_list.end())
      {
	i_next++;
      }
    if(i!=section_list.begin())
      {
	i_prev--;
      }
    if(i==section_list.begin())
      {
	T x1=i->x_l();
	T x2=i->x_l();
	T x3=i_next->x_l();
	T x4=i_next->x_u();

	T y1=i->y_l();
	T y2=std::numeric_limits<T>::infinity();
	T y3=i_next->y_l();
	T y4=i_next->y_u();

	if(x3==x4)
	  {
	    throw ill_conditioned_distribution();
	  }
	  
	
#ifdef DEBUG
	assert(!std::isnan(y1));
	assert(!std::isnan(y2));
	assert(!std::isnan(y3));
	assert(!std::isnan(y4));

#endif
	auto p=solve_intersection(*i,
				  std::make_pair(x1,y1),
				  std::make_pair(x2,y2),
				  std::make_pair(x3,y3),
				  std::make_pair(x4,y4));
	i->set_x_i(p.first);
	i->set_y_i(p.second);
      }
    else if(i_next==section_list.end())
      {
	T x1=i_prev->x_l();
	T x2=i_prev->x_u();
	T x3=i->x_u();
	T x4=i->x_u();

	T y1=i_prev->y_l();
	T y2=i_prev->y_u();
	T y3=i->y_u();
	T y4=std::numeric_limits<T>::infinity();

	if(x1==x2)
	  {
	    throw ill_conditioned_distribution();
	  }

#ifdef DEBUG
	assert(!std::isnan(y1)&&
	       !std::isnan(y2)&&
	       !std::isnan(y3)&&
	       !std::isnan(y4));
#endif
	auto p=solve_intersection(*i,
				  std::make_pair(x1,y1),
				  std::make_pair(x2,y2),
				  std::make_pair(x3,y3),
				  std::make_pair(x4,y4));
	i->set_x_i(p.first);
	i->set_y_i(p.second);
      }
    else
      {
	T x1=i_prev->x_l();
	T x2=i_prev->x_u();
	T x3=i_next->x_l();
	T x4=i_next->x_u();

	T y1=i_prev->y_l();
	T y2=i_prev->y_u();
	T y3=i_next->y_l();
	T y4=i_next->y_u();

	if(x1==x2||x3==x4)
	  {
	    throw ill_conditioned_distribution();
	  }

#ifdef DEBUG
	assert(!std::isnan(y1)&&
	       !std::isnan(y2)&&
	       !std::isnan(y3)&&
	       !std::isnan(y4));
#endif
	auto p=solve_intersection(*i,
				  std::make_pair(x1,y1),
				  std::make_pair(x2,y2),
				  std::make_pair(x3,y3),
				  std::make_pair(x4,y4));
	i->set_x_i(p.first);
	i->set_y_i(p.second);
      }
    if(i->x_i() < i->x_l())
      {
	i->set_x_i( i->x_l() );
      }
    if(i->x_i() > i->x_u())
      {
	i->set_x_i( i->x_u() );
      }
#ifdef DEBUG
    assert(i->x_i() >= i->x_l());
    assert(i->x_i() <= i->x_u());
    assert(!std::isinf(i->y_i())||i->y_i()<0);
#endif
  }

  template <typename T>
  T calc_scale(const std::list<section<T> >& section_list)
  {
    T scale=-std::numeric_limits<T>::infinity();
    for(auto& i:section_list)
      {
	scale=std::max(scale,std::max(i.y_l(),i.y_u()));
      }
    if(std::isinf(scale))
      {
	if(scale<0)
	  {
	    throw ill_conditioned_distribution("maybe all values are zero 7");
	  }
	else
	  {
	    throw ill_conditioned_distribution("may have inf points 8");
	  }
      }
    
    return scale;
  }

  template <typename T>
  void update_scale(std::list<section<T> >& section_list,T& scale)
  {
    T new_scale=calc_scale(section_list);
    for(auto i=section_list.begin();i!=section_list.end();++i)
      {
	i->set_y_l(i->y_l()-new_scale);
	i->set_y_u(i->y_u()-new_scale);
      }
    
    scale+=new_scale;
    for(auto i=section_list.begin();i!=section_list.end();++i)
      {
	calc_intersection(section_list,i);
	calc_cum_int_exp_y(section_list,i);
#ifdef DEBUG
	assert(i->x_l()<=i->x_i()&&
	       i->x_i()<=i->x_u());
#endif
      }
#ifdef DEBUG
    assert(!std::isnan(section_list.back().cum_int_exp_y_u()));
#endif
  }


  template <typename T,typename TD>
  void init(const TD& pd,const std::pair<T,T>& xrange, const std::vector<T>& init_x1, std::list<section<T> >& section_list,T& scale)
  {
    if(init_x1.size()<3)
      {
	throw too_few_init_points();
      }

    for(size_t i=1;i<init_x1.size();++i)
      {
	if(init_x1[i]<=init_x1[i-1])
	  {
	    throw data_not_in_order();
	  }
      }
    std::vector<T> init_x{xrange.first};
    for(auto& x:init_x1)
      {
	if(std::abs(x-init_x.back())<std::numeric_limits<T>::epsilon()*10||
	   x<=init_x.back())
	  {
	    continue;
	  }
	init_x.push_back(x);
      }
    //init_x.push_back(xrange.second-std::numeric_limits<T>::epsilon());
    init_x.push_back(xrange.second);
    if(init_x.size()<5)
      {
	too_few_init_points e;
	e.attach_message("#8");
	throw e;
      }
    
    for(size_t i=0;i!=init_x.size()-1;++i)
      {
	section<T> s;
	
	s.set_x_l(init_x[i]);
	s.set_x_u(init_x[i+1]);
	s.set_y_l(eval_log(pd,init_x[i],(T)0));
	s.set_y_u(eval_log(pd,init_x[i+1],(T)0));
	if(!(xrange.first<=s.x_l()&&s.x_l()<=s.x_u()&&s.x_u()<=xrange.second))
	  {
	    data_not_in_order e;
	    e.attach_message("#9");
	    throw e;
	  }

	if(std::isnan(s.y_l())||std::isnan(s.y_u()))
	  {
#ifdef DEBUG
	    std::cerr<<s.x_l()<<" "<<s.x_u()<<std::endl;
	    std::cerr<<s.y_l()<<" "<<s.y_u()<<std::endl;
    	    assert(0);
#endif
	    nan_or_inf e;
	    e.attach_message("#10");
	    throw e;
	  }
	
	section_list.push_back(s);
      }
    
    scale=calc_scale(section_list);
    
    //scale=0;
    for(auto i=section_list.begin();i!=section_list.end();++i)
      {
	i->set_y_l(i->y_l()-scale);
	i->set_y_u(i->y_u()-scale);
      }

    
    for(auto i=section_list.begin();i!=section_list.end();++i)
      {
	calc_intersection(section_list,i);
	calc_cum_int_exp_y(section_list,i);
#ifdef DEBUG
	assert(i->x_l()<=i->x_i()&&
	       i->x_i()<=i->x_u());
#endif
      }

#ifdef DEBUG
    assert(!std::isnan(section_list.back().cum_int_exp_y_u()));
#endif

    for(size_t i=0;;++i)
      {
	assert(i<1000);
	bool has_inf=false;
	auto iter=section_list.begin();
	
	for(;iter!=section_list.end();++iter)
	  {
	    if(std::isinf(iter->cum_int_exp_y_u()))
	      {
		has_inf=true;
		T x=(iter->x_l()+iter->x_u())/2;

		insert_point(pd,section_list,x,scale);
		update_scale(section_list,scale);
		break;
	      }
	  }
	
	if(!has_inf)
	  {
	    break;
	  }
      }
    

    if(std::isinf(section_list.back().cum_int_exp_y_u()))
      {
#ifdef DEBUG
	std::cerr<<"initial points:"<<std::endl;
	for(auto & i:section_list)
	  {
	    std::cerr<<i.x_l()<<" "<<i.y_l()<<"|"<<i.x_i()<<" "<<i.y_i()<<"|"<<i.x_u()<<" "<<i.y_u()<<std::endl;
	  }
#endif
	more_init_points_needed e;
	e.attach_message("#11");
	throw e;
      }
  }

  enum class insertion_result
  {
    SUCCEEDED,SEARCH_FAILED,POINT_OVERLAPPED
  };

  template <typename T,typename TD>
  insertion_result insert_point(const TD& pd,std::list<section<T> >& section_list,const T& x,T scale)
  {
#if 0
    auto iter=section_list.begin();
    
    for(;iter!=section_list.end();++iter)
      {
	if(iter->encloses(x)&&x!=iter->x_l()&&x!=iter->x_u())
	  {
	    break;
	  }
      }

#else
    auto iter=lower_bound(section_list.begin(),section_list.end(),x,[](const section<T>& x,const T& y)
		       {
			 if( (x.x_l()-y)*(x.x_u()-y)<=0  )
			   {
			     return false;
			   }
			 return x.x_l()<y;
		       });
#endif    
      
    if(iter==section_list.end())
      {
	return insertion_result::SEARCH_FAILED;
      }
#if 0
    if(std::abs(iter->x_l()-x)<std::numeric_limits<T>::epsilon()*10||std::abs(iter->x_u()-x)<std::numeric_limits<T>::epsilon()*10)
      {
	return insertion_result::POINT_OVERLAPPED;
      }

    
    if(std::abs(x-iter->x_l())<std::numeric_limits<T>::epsilon()*10||std::abs(x-iter->x_u())<std::numeric_limits<T>::epsilon()*10)
      {
	return insertion_result::POINT_OVERLAPPED;
      }
#endif
    section<T> s;
    s.set_x_l(iter->x_l());
    s.set_x_u(x);
    s.set_y_l(eval_log(pd,s.x_l(),scale));
    s.set_y_u(eval_log(pd,s.x_u(),scale));

    iter->set_x_l(x);
    iter->set_y_l(eval_log(pd,iter->x_l(),scale));

#ifdef DEBUG
    assert(!std::isnan(iter->y_l()));
    assert(!std::isnan(iter->y_u()));
#endif
    section_list.insert(iter,s);
    return insertion_result::SUCCEEDED;
  }



  template <typename T,typename TD>
  typename std::list<section<T> >::const_iterator search_point(const std::list<section<T> >& section_list,T p,
							       const TD& pd,const T scale)
  {
    /*
    section<T> v;
    v.set_cum_int_exp_y_u(p*section_list.back().cum_int_exp_y_u());
    
    auto r=equal_range(section_list.begin(),section_list.end(),v,[](const section<T>& x,const section<T>& y){return x.cum_int_exp_y_u()<y.cum_int_exp_y_u();});

    if(r.first==section_list.end()||r.second==section_list.end())
      {
	std::cerr<<(p*section_list.back().cum_int_exp_y_u())<<std::endl;
	std::cerr<<section_list.back().cum_int_exp_y_u()<<std::endl;
	throw search_failed();
      }
    while(1)
      {
	if(r.first==r.second)
	  {
	    break;
	  }
	if(++r.first==r.second)
	  {
	    break;
	  }
	if(r.first==--r.second)
	  {
	    break;
	  }
      }
    
    assert(r.first!=section_list.end());
    
    if(!(r.first->x_l()<=r.first->x_i()&&
	 r.first->x_i()<=r.first->x_u()))
      {
	std::cerr<<r.first->x_l()<<" "<<r.first->x_i()<<" "<<r.first->x_u()<<std::endl;
	assert(0);
      }
    */
    if(section_list.back().cum_int_exp_y_u()==0)
      {
	throw ill_conditioned_distribution("maybe all values are zero #1059");
      }
    
    T x=p*section_list.back().cum_int_exp_y_u();

    for(auto i=section_list.begin();i!=section_list.end();++i)
      {
	const T xx=i->cum_int_exp_y_u();
#ifdef DEBUG
	assert(!std::isnan(xx));
#endif

	if((p<1&&x<xx)||(p==1&&x<=xx))
	  {
	    return i;
	  }
      }
    
#ifdef DEBUG
    std::ofstream ofs_log("log.txt");
    ofs_log<<std::setprecision(20);
    ofs_log<<section_list;
    ofs_log.close();

    ofs_log.open("dump.txt");
    ofs_log<<std::setprecision(20);
    for(auto& i:section_list)
      {
	ofs_log<<i.x_l()<<" "<<i.y_l()<<std::endl;
	ofs_log<<i.x_i()<<" "<<i.y_i()<<std::endl;
	ofs_log<<i.x_u()<<" "<<i.y_u()<<std::endl;
      }
    ofs_log<<"no no no"<<std::endl;
    for(auto& i:section_list)
      {
	ofs_log<<i.x_l()<<" "<<eval_log(pd,i.x_l(),scale)<<std::endl;
	ofs_log<<i.x_i()<<" "<<eval_log(pd,i.x_i(),scale)<<std::endl;
	ofs_log<<i.x_u()<<" "<<eval_log(pd,i.x_u(),scale)<<std::endl;
      }
    ofs_log.close();
    assert(0);
#endif

    search_failed e;
    e.attach_message("#12");
    throw e;
    return section_list.end();
  }

  template <typename T,typename TD>
  void check_range(const TD& pd,std::list<section<T> >& section_list,T& scale)
  {
    for(size_t i=0;;++i)
      {
	bool has_inf=false;
	auto iter=section_list.begin();
	for(;iter!=section_list.end();++iter)
	  {
	    if(std::isinf(iter->cum_int_exp_y_l()))
	      {
		has_inf=true;
		T x=(iter->x_l()+iter->x_i())/2;
		if(insert_point(pd,section_list,x,scale)==insertion_result::SUCCEEDED)
		  {
		    update_scale(section_list,scale);
		    break;
		  }
		else
		  {
		    throw ill_conditioned_distribution();
		  }
	      }
	    else if(std::isinf(iter->cum_int_exp_y_u()))
	      {
		has_inf=true;
		T x=(iter->x_i()+iter->x_u())/2;
		
		if(insert_point(pd,section_list,x,scale)==insertion_result::SUCCEEDED)
		  {
		    update_scale(section_list,scale);
		    break;
		  }
		else
		  {
		    throw ill_conditioned_distribution();
		  }
	      }
	  }
	
	if(!has_inf)
	  {
	    break;
	  }
      }
  }

  
  template <typename T,typename T_urand,typename TD>
  T sample(const std::list<section<T> >& section_list,T_urand& rnd,const TD& pd,T scale)
  {
    T result=0;

#if 0
    static int n=0;
    
    std::string fname("cum_");
    fname+=std::to_string(n);
    fname+=".qdp";
    n++;
    std::ofstream ofs(fname.c_str());
    ofs<<std::setprecision(20);
    for(auto& i:section_list)
      {
	ofs<<i.x_i()<<" "<<i.cum_int_exp_y_l()<<std::endl;
	ofs<<i.x_u()<<" "<<i.cum_int_exp_y_u()<<std::endl;
      }
#endif
    T p=0;
    for(size_t i=0;i<10;++i)
      {
	p=rnd();
	if(p>=0&&p<1)
	  {
	    break;
	  }
	else
	  {
#ifdef DEBUG
	    if(i>5)
	      {
		assert(0);
	      }
#endif
	  }
      }
    //std::cerr<<"p="<<p<<" ";
    auto iter=search_point(section_list,p,pd,scale);
    T y=section_list.back().cum_int_exp_y_u()*p;
    
#ifdef DEBUG
    assert(!std::isnan(p));
    assert(!std::isnan(section_list.back().cum_int_exp_y_u()));
    
    if(std::isnan(y))
      {
	std::cerr<<"p="<<p<<std::endl;
	std::cerr<<"cum_y="<<section_list.back().cum_int_exp_y_u()<<std::endl;
      }
    
    assert(!std::isnan(y));
#endif
    
    T x1,x2,y1,y2;
    T ybase;
    
    if(y>=iter->cum_int_exp_y_l()&&
       y<=iter->cum_int_exp_y_u())
      {
	ybase = iter->cum_int_exp_y_l();
	x1 = iter -> x_i();
	y1 = iter -> y_i();
	
	x2 = iter -> x_u();
	y2 = iter -> y_u();
	
#ifdef DEBUG
	assert(x2>=x1);
#endif
      }
    else if(y<iter->cum_int_exp_y_l())
      {
	auto iter1=iter;
	
	if(iter1!=section_list.begin())
	  {
	    --iter1;
	    ybase = iter1->cum_int_exp_y_u();
	  }
	else
	  {
	    ybase = 0;
	  }
	
	x1 = iter -> x_l();
	y1 = iter -> y_l();
	
	x2 = iter -> x_i();	   
	y2 = iter -> y_i();
#ifdef DEBUG	    
	assert(x2>=x1);
#endif
      }
    else
      {
	assert(0);
      }
    
    if ( y == ybase )
      {
	result = x1;
      }
    else
      {
#ifdef DEBUG
	assert(!std::isnan(y));
	assert(!std::isinf(y));
	assert(!std::isnan(ybase));
	assert(!std::isinf(ybase));
	assert(x2>=x1);
#endif
	result=inv_int_exp_y(y-ybase,std::make_pair(x1,y1),std::make_pair(x2,y2));
      }
    if(std::isnan(result))
      {
#ifdef DEBUG	  
	assert(0);
#endif
	nan_or_inf e;
	e.attach_message("#13");
	throw e;
      }
    //std::cerr<<result<<std::endl;
    return result;
  }
  
  template <typename T,typename TD,typename T_urand>
  T arms(const TD& pd,const std::pair<T,T>& xrange,
	 const std::vector<T>& init_x,T xcur,size_t n,T_urand& rnd,size_t& xmchange_count)
  {
    if(xrange.second<xrange.first+std::numeric_limits<T>::epsilon()*std::abs(xrange.first))
      {
	throw too_small_var_range();
      }
    if(xcur<xrange.first||xcur>xrange.second)
      {
#ifdef DEBUG
	std::cerr<<"xcur="<<xcur<<" ("<<xrange.first<<" , "<<xrange.second<<")"<<std::endl;
#endif
	var_out_of_range e;
	e.attach_message("#14");
	throw e;
      }

    std::list<section<T> > section_list;
    
    T scale=0;

    init(pd,xrange,init_x,section_list,scale);
    
    T xm=-1;
    //bool xmchanged=false;
    //size_t xmchange_count=0;
    //assert(xcur>=xrange.first&&xcur<=xrange.second);

    
    for(size_t i=0,cnt=0;i<n;)
      {
	T x=0;
	for(size_t j=0;;++j)
	  {
	    
	    x=sample(section_list,rnd,pd,scale);
	    
	    if(j>10000)
	      {
		if(x<=xrange.first)
		  {
		    x=xrange.first+std::numeric_limits<T>::epsilon()*std::abs(xrange.first);
		  }
		else if(x>=xrange.second)
		  {
		    x=xrange.second-std::numeric_limits<T>::epsilon()*std::abs(xrange.second);
		  }
		break;
	      }
	    
	    if(x<=xrange.first)
	      {
		T x1=x;
		for(auto & s:section_list)
		  {
		    if (x<s.x_i())
		      {
			x1=s.x_i();
			break;
		      }
		    if(x<s.x_u())
		      {
			x1=s.x_u();
			break;
		      }
		  }
		x1=(x1+x)/2;
		insert_point(pd,section_list,x1,scale);
		update_scale(section_list,scale);
		if(std::isinf(section_list.back().cum_int_exp_y_u()))
		  {
		    check_range(pd,section_list,scale);
		  }
		continue;
	      }
	    if(x>=xrange.second)
	      {
		T x1=x;
		for(auto p=section_list.rbegin();p!=section_list.rend();++p)
		  {
		    if(p->x_i()<x)
		      {
			x1=p->x_i();
			break;
		      }
		    if(p->x_l()<x)
		      {
			x1=p->x_l();
			break;
		      }
		  }
		x1=(x1+x)/2;
		insert_point(pd,section_list,x1,scale);
		update_scale(section_list,scale);
		
		if(std::isinf(section_list.back().cum_int_exp_y_u()))
		  {
		    check_range(pd,section_list,scale);
		  }
		continue;
	      }
	    
	    break;
	  }
	
	T u=rnd();
	T xa=0;
	if(std::log(u)+eval(x,section_list)>eval_log(pd,x,scale))
	  {
	    insert_point(pd,section_list,x,scale);
	    update_scale(section_list,scale);
	    if(std::isinf(section_list.back().cum_int_exp_y_u()))
	      {
		check_range(pd,section_list,scale);
	      }
#ifdef DEBUG
	    assert(!std::isinf(section_list.back().cum_int_exp_y_u()));
	    assert(!std::isnan(section_list.back().cum_int_exp_y_u()));
#endif
	    if(cnt++>100*i)
	      {
		
#if 0
		std::ofstream ofs("dump.qdp");
		ofs<<std::setprecision(20);
		for(auto& i:section_list)
		  {
		    ofs<<i.x_l()<<" "<<i.y_l()<<std::endl;
		    ofs<<i.x_i()<<" "<<i.y_i()<<std::endl;
		    ofs<<i.x_u()<<" "<<i.y_u()<<std::endl;	
		  }
		
		ofs<<"no no no"<<std::endl;
		
		for(T x=0.001;x<100;x+=.01)
		  {
		    ofs<<x<<" "<<eval_log(pd,x,scale)<<std::endl;
		  }
		ofs.close();
		ofs.open("cum.qdp");
		for(auto& i:section_list)
		  {
		    ofs<<i.x_i()<<" "<<i.cum_int_exp_y_l()<<std::endl;
		    ofs<<i.x_u()<<" "<<i.cum_int_exp_y_u()<<std::endl;
		  }
		//exit(0);
#endif
		if(cnt>100)
		  {
		    throw too_many_rejections();
		  }
	      }
	    continue;
	  }
	else
	  {
	    xa=x;
	  }
	u=rnd();
	
	if(std::log(u)>std::min(static_cast<T>(0),eval_log(pd,xa,scale)-eval_log(pd,xcur,scale)+std::min(eval_log(pd,xcur,scale),eval(xcur,section_list))-std::min(eval_log(pd,xa,scale),eval(xa,section_list))))
	  {
	    xm=xcur;
	    ++i;
	  }
	else
	  {
	    ++xmchange_count;
	    xm=xa;
	    ++i;
	  }
	xcur=xm;
      }

    
    return xm;
  }


  template <typename T,typename TD,typename T_urand>
  T arms(const TD& pd,const std::pair<T,T>& xrange,
	 T xcur,size_t n,T_urand& rnd,size_t& xmchange_count)
  {
    std::vector<T> init_x(5);
    T xl=xrange.first,xr=xrange.second;
    for(size_t n=0;n<init_x.size();++n)
	{
	  init_x[n]= xl+(xr-xl)/(init_x.size()+1)*(n+1);
	}
    return arms(pd,xrange,init_x,xcur,n,rnd,xmchange_count);
  }
}


#endif
