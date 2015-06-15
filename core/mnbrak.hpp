#ifndef MCMC_MNBRAK_HPP
#define MCMC_MNBRAK_HPP 

//#include "optimizer.hpp"
#include "bas_util.hpp"
namespace mcmc_utilities
{

  
  template <typename T_dist>
  void mnbrak(typename T_dist::T_var& ax,typename T_dist::T_var& bx,typename T_dist::T_var& cx,typename T_dist::T_p& fa,typename T_dist::T_p& fb,typename T_dist::T_p& fc,const dist_adapter<typename T_dist::T_p,typename T_dist::T_var>& func)
  {
    const typename T_dist::T_var GOLD=1.618034;
    const typename T_dist::T_var GLIMIT=100;
    const typename T_dist::T_var TINY=std::numeric_limits<typename T_dist::T_var>::epsilon();
    typename T_dist::T_var ulim,u,r,q,fu;
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
	u=bx-typename T_dist::T_var((bx-cx)*q-(bx-ax)*r)/
	  typename T_dist::T_var(typename T_dist::T_var(2.)*sign(typename T_dist::T_var(std::max(typename T_dist::T_var(std::abs(typename T_dist::T_var(q-r))),typename T_dist::T_var(TINY))),typename T_dist::T_var(q-r)));
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
		shft3(bx,cx,u,typename T_dist::T_var(cx+GOLD*(cx-bx)));
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
