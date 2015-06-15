#ifndef MCMC_BRENT_HPP
#define MCMC_BRENT_HPP

//#include <iostream>
#include "bas_util.hpp"
//#include "optimizer.hpp"
namespace mcmc_utilities
{
  template<typename T_dist>
  typename T_dist::T_p brent(typename T_dist::T_var ax,typename T_dist::T_var bx,typename T_dist::T_var cx,const T_dist& f,typename T_dist::T_var tol,typename T_dist::T_var& xmin)
  {
    const int ITMAX=100;
    const typename T_dist::T_var CGOLD=0.3819660;
    const typename T_dist::T_var ZEPS=std::numeric_limits<typename T_dist::T_p>::epsilon()*1.e-3;
    
    int iter;
    typename T_dist::T_var a=0,b=0,d(0),etemp=0,p=0,q=0
      ,r=0,tol1=0,tol2=0,u=0,v=0,w=0,x=0,xm=0;
    typename T_dist::T_p fu=0,fv=0,fw=0,fx=0;
    typename T_dist::T_p e=0.;
    a=(ax<cx?ax:cx);
    b=(ax>cx?ax:cx);
    x=w=v=bx;
    fw=fv=fx=f(x);
    for(iter=0;iter<ITMAX;++iter)
      {
	xm=.5*(a+b);
	tol2=2.*(tol1=tol*std::abs(x)+ZEPS);
	if(std::abs(typename T_dist::T_var(x-xm))<=(tol2-.5*(b-a)))
	  {
	    xmin=x;
	    return fx;
	  }
	if(std::abs(e)>tol1)
	  {
	    r=(x-w)*(fx-fv);
	    q=(x-v)*(fx-fw);
	    p=(x-v)*q-(x-w)*r;
	    q=2.*(q-r);
	    if(q>0.)
	      {
		p=-p;
	      }
	    q=std::abs(q);
	    etemp=e;
	    e=d;
	    if(std::abs(p)>=std::abs(typename T_dist::T_p(typename T_dist::T_p(.5)*p*etemp))||p<=q*(a-x)||p>=q*(b-x))
	      {
		d=CGOLD*(e=(x>=xm?a-x:b-x));
	      }
	    else
	      {
		d=p/q;
		u=x+d;
		if(u-a<tol2||b-u<tol2)
		  {
		    d=sign(tol1,typename T_dist::T_var(xm-x));
		  }
	      }
	    
	  }
	else
	  {
	    d=CGOLD*(e=(x>=xm?a-x:b-x));
	  }
	u=(std::abs(d)>=tol1?x+d:x+sign(tol1,d));
	fu=f(u);
	if(fu<=fx)
	  {
	    if(u>=x)
	      {
		a=x;
	      }
	    else
	      {
		b=x;
	      }
	    shft3(v,w,x,u);
	    shft3(fv,fw,fx,fu);
	  }
	else
	  {
	    if(u<x)
	      {
		a=u;
	      }
	    else
	      {
		b=u;
	      }
	    if(fu<=fw||w==x)
	      {
		v=w;
		w=u;
		fv=fw;
		fw=fu;
	      }
	    else if(fu<=fv||v==x||v==w)
	      {
		v=u;
		fv=fu;
	      }
	  }
      }
    //std::cerr<<"Too many iterations in brent"<<std::endl;
    xmin=x;
    return fx;
    
  }
}

#endif
