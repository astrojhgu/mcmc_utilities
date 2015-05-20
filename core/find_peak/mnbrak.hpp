#ifndef MCMC_MNBRAK_HPP
#define MCMC_MNBRAK_HPP 

//#include "optimizer.hpp"
#include "bas_util.hpp"
namespace mcmc_utilities
{

  
  template <typename T_p,typename T_var>
  void mnbrak(T_var& ax,T_var& bx,T_var& cx,T_p& fa,T_p& fb,T_p& fc,const cfunc<T_p,T_var>& func)
  {
    const T_var GOLD=1.618034;
    const T_var GLIMIT=100;
    const T_var TINY=std::numeric_limits<T_var>::epsilon();
    T_var ulim,u,r,q,fu;
    fa=func(ax);
    fb=func(bx);
    
    if(fb>fa)
      {
	//shft(dum,ax,bx,dum);
	//shft(dum,fb,fa,dum);
	std::swap(ax,bx);
	std::swap(fa,fb);
      }

    cx=bx+GOLD*(bx-ax);
    fc=func(cx);
    while(fb>fc)
      {
	r=(bx-ax)*(fb-fc);
	q=(bx-cx)*(fb-fa);
	u=bx-T_var((bx-cx)*q-(bx-ax)*r)/
	  T_var(T_var(2.)*sign(T_var(tmax(T_var(tabs(T_var(q-r))),T_var(TINY))),T_var(q-r)));
	ulim=bx+GLIMIT*(cx-bx);
	if((bx-u)*(u-cx)>0.)
	  {
	    fu=func(u);
	    if(fu<fc)
	      {
		ax=bx;
		bx=u;
		fa=fb;
		fb=fu;
		return;
	      }
	    else if(fu>fb)
	      {
		cx=u;
		fc=fu;
		return;
	      }
	    u=cx+GOLD*(cx-bx);
	    fu=func(u);
	  }
	else if((cx-u)*(u-ulim)>0.)
	  {
	    fu=func(u);
	    if(fu<fc)
	      {
		shft3(bx,cx,u,T_var(cx+GOLD*(cx-bx)));
		shft3(fb,fc,fu,func(u));
	      }
	  }
	else if((u-ulim)*(ulim-cx)>=0)
	  {
	    u=ulim;
	    fu=func(u);
	  }
	else
	  {
	    u=cx+GOLD*(cx-bx);
	    fu=func(u);
	  }
	shft3(ax,bx,cx,u);
	shft3(fa,fb,fc,fu);
      }
  }
}

#endif
