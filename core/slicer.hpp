#ifndef MYSLICER_HPP
#define MYSLICER_HPP
#include "distribution.hpp"
#include "mcmc_exception.hpp"
#include "base_urand.hpp"
#include <algorithm>
//#define DEBUG
#ifdef DEBUG
#include <fstream>
#endif
#include <cassert>
#include <list>
#include <utility>
#include <limits>
#include <cmath>

namespace mcmc_utilities
{

  template <typename T>
  class slice_sampler
  {
  private:
    bool adapt;
    T width;
    size_t nmax;
    size_t niter;
    T sumdiff;
    const size_t min_adapt;
    const probability_density_1d<T>& pd;
    
  public:
    slice_sampler(const probability_density_1d<T>& _pd, const T& _width,size_t _nmax)
      :adapt(true),width(_width),nmax(_nmax),niter(0),sumdiff(0),min_adapt(50),pd(_pd)
    {
    }
    
  private:
    T exponential(const base_urand<T>& rnd)
    {
      T result=0;
      do
	{
	  result=-std::log(1-rnd());
	}
      while(std::isinf(result));
      return result;
    }
    
    bool accept(T xold, T xnew, T z, T L, T R,
		T lower, T upper)
    {
      //Acceptance step for doubling update method
      
      bool d = false;
      while ((R - L) > 1.1 * width)
	{
	  T M = (L + R)/2;
	  if ((xold < M && xnew >= M) || (xold >= M && xnew < M))
	    d = true;
	  if (xnew < M)
	    {
	      R = M;
	    }
	  else
	    {
	      L = M;
	    }
	  if (d)
	    {
	      bool right_ok = true;
	      if (R <= upper)
		{
		  right_ok = pd.eval_log(R) < z;
		}
	      bool left_ok = true;
	      if (L >= lower)
		{
		  left_ok = pd.eval_log(L) < z;
		}
	      if (left_ok && right_ok)
		{
		  return false;
		}
	    }
	}
      return true;
    }
    
  public:
    T sample(T& xcur,const base_urand<T>& rng)
    {
      T g0 = pd.eval_log(xcur);
      if (std::isnan(g0)||std::isinf(g0))
	{
	  throw nan_or_inf();
	}
      
      // Generate auxiliary variable
      T z = g0 - exponential(rng);
      
      // Generate random interval of width "_width" about current value
      T xold = xcur;
      T L = xold - rng() * width; 
      T R = L + width;
      
      T lower = 0;
      T upper = 0;
      
      auto px=pd.var_range();
      lower=px.first;
      upper=px.second;
      
      // Doubling 
      bool left_ok = false, right_ok = false;
      for (size_t i = 0; i < nmax; ++i) {
	if (rng() < 0.5)
	  {
	    if (L >= lower)
	      {
		L = 2*L - R;
		if (L < lower)
		  {
		    left_ok = true;
		  }
		else
		  {
		    //setValue(L);
		    xcur=L;
		    left_ok = pd.eval_log(xcur) < z;
		  }
	      }
	    else
	      {
		left_ok = true;
	      }
	  }
	else
	  {
	    if (R <= upper)
	      {
		R = 2*R - L;
		if (R > upper)
		  {
		    right_ok = true;
		  }
		else
		  {
		    //setValue(R);
		    xcur=R;
		    right_ok = pd.eval_log(xcur) < z;
		  }
	      }
	    else
	      {
		right_ok = true;
	      }
	  }
	if (left_ok && right_ok)
	  {
	    break;
	  }
      }
      
      // Keep sampling from the interval until acceptance (the loop is
      // guaranteed to terminate).
      T Lbar = L, Rbar = R;
      T xnew;
      for(;;)
	{
	  xnew =  Lbar + rng() * (Rbar - Lbar);
	  if (xnew >= lower && xnew <= upper)
	    {
	      xcur=xnew;
	      T g = pd.eval_log(xnew);
	      if (g >= z && accept(xold, xnew, z, L, R, lower, upper))
		{
		  //setValue(xnew);
		  xcur=xnew;
		  break;
		}
	    }
	  // shrink the interval
	  if (xnew <= xold)
	    {
	      Lbar = xnew;
	    }
	  else
	    {
	      Rbar = xnew;
	    }
	}

      if (adapt)
	{
	  sumdiff += niter * std::abs(xnew - xold);
	  ++niter;
	  if (niter > min_adapt)
	    {
	      width = 2 * sumdiff / niter / (niter - 1);  
	    }
	}
      
      return xcur;
    }
  };
}

#endif
