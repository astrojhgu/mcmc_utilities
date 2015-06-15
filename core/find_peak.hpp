/**
   \file linmin.hpp
   \brief linear search
   \author Junhua Gu
 */


#ifndef MCMC_LINMIN_HPP
#define MCMC_LINMIN_HPP

#include "mnbrak.hpp"
#include "brent.hpp"
#include <cmath>
#include "mcmc_traits.hpp"
//#include <iostream>
using namespace std;
#include <limits>
namespace mcmc_utilities
{
  template<typename T_p,typename T_var>
  T_var find_peak(const probability_density_1d<T_p,T_var>& dist)
  {
    //func_adaptor<T_p,T_var> fadpt(p,xi,func);
    typedef dist_adapter<T_p,T_var> T_dist;
    T_dist func;
    
    func.ppd=&dist;
    dist.var_range(func.xl,func.xr);

    int j=0;
    const T_p TOL=std::sqrt(std::numeric_limits<T_p>::epsilon());
    T_p xx=0,peak_x=0,fx=0,fb=0,fa=0,bx=0,ax=0;



    ax=func.xl;
    
    bx=func.xr;
    xx=(ax+bx)/2;

    mnbrak<T_dist>(ax,xx,bx,fa,fx,fb,func);
    //cout<<xx<<endl;
    brent<T_dist>(ax,xx,bx,func,TOL,peak_x);

    return peak_x;
  }
}


#endif
