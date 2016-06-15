#ifndef MYSLICER_HPP
#define MYSLICER_HPP
#include "distribution.hpp"
#include "error_handler.hpp"
#include "base_urand.hpp"
#include <algorithm>
#include <limits>
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

  template <typename T, typename TD>
  class slice_sampler
  {
  private:
    bool adapt;
    T width;
    size_t nmax;
    size_t nmin;
    size_t niter;
    T sumdiff;
    const size_t min_adapt;
    const TD& pd;
    const std::pair<T,T>& xrange;
    
  public:
    slice_sampler(const TD& _pd, const std::pair<T,T>& _xrange, const T& _width,size_t _nmax,size_t _nmin)
      :adapt(true),width(_width),nmax(_nmax),nmin(_nmin),niter(0),sumdiff(0),min_adapt(50),pd(_pd),xrange(_xrange)
    {
    }

    slice_sampler(const slice_sampler<T,TD>&)=delete;
    slice_sampler<T,TD>& operator=(const slice_sampler<T,TD>&)=delete;
    
  private:
    T exponential(base_urand<T>& rnd)const
    {
      T result=0;
      do
	{
	  result=-std::log(1-rnd());
	}
      while(std::isinf(result));
      return result;
    }
    
    bool accept(const T& xold,const T& xnew,const T& z,T L, T R,
		const T& lower, const T& upper)const
    {
      bool d = false;
      while ((R - L) > static_cast<T>(1.1) * width)
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
		  right_ok = pd(R) < z;
		}
	      bool left_ok = true;
	      if (L >= lower)
		{
		  left_ok = pd(L) < z;
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
    T sample_double(T& xcur,base_urand<T>& rng)
    {
      for(size_t i=0;i<nmin;++i)
	{
	  sample1_double(xcur,rng);
	}
      return xcur;
    }
    
    T sample1_double(T& xcur,base_urand<T>& rng)
    {
      T lower = 0;
      T upper = 0;
      
      //auto xrange=pd.var_range();
      lower=xrange.first;
      upper=xrange.second;

      if(xcur<xrange.first||xcur>xrange.second)
      {
	std::cerr<<"xcur="<<xcur<<" ("<<xrange.first<<" , "<<xrange.second<<")"<<std::endl;
	throw var_out_of_range();
      }

      
      T g0 = pd(xcur);
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
      
      
      // Doubling 
      bool left_ok = false, right_ok = false;
      for (size_t i = 0; i < nmax; ++i) {
	if (rng() < static_cast<T>(0.5))
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
		    xcur=L;
		    left_ok = pd(xcur) < z;
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
		    xcur=R;
		    right_ok = pd(xcur) < z;
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
      
      T Lbar = L, Rbar = R;
      T xnew;
      for(;;)
	{
	  xnew =  Lbar + rng() * (Rbar - Lbar);
	  if (xnew >= lower && xnew <= upper)
	    {
	      xcur=xnew;
	      T g = pd(xnew);
	      if (g >= z && accept(xold, xnew, z, L, R, lower, upper))
		{
		  xcur=xnew;
		  break;
		}
	    }

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

    T sample_step(T& xcur,base_urand<T>& rng)
    {
      for(size_t i=0;i<nmin;++i)
	{
	  sample1_step(xcur,rng);
	}
      return xcur;
    }

    T sample1_step(T& xcur,base_urand<T>& rng)
    {
      T lower = 0;
      T upper = 0;
      //auto xrange=pd.var_range();
      lower=xrange.first;
      upper=xrange.second;
      if(xcur<xrange.first||xcur>xrange.second)
	{
	  std::cerr<<"xcur="<<xcur<<" ("<<xrange.first<<" , "<<xrange.second<<")"<<std::endl;
	  throw var_out_of_range();
	}

      T g0 = pd(xcur);
      if (std::isinf(g0))
	{
	  std::cerr<<"xcur="<<xcur<<std::endl;
	  throw nan_or_inf();
	}

      // Generate auxiliary variable
      T z = g0 - exponential(rng);
      
      // Generate random interval of width "_width" about current value
      T xold = xcur;
      T L = xold - rng() * width; 
      T R = L + width;
      
      size_t j = static_cast<size_t>(rng() * nmax);
      size_t k = nmax - 1 - j;
      
      
      if (L < lower)
	{
	  L = lower;
	}
      else
	{
	  //setValue(L);
	  xcur=L;
	  while (j-- > 0 && pd(xcur) > z)
	    {
	      L -= width;
	      if (L < lower)
		{
		  L = lower;
		  break;
		}
	      //setValue(L);
	      xcur=L;
	    }
	}

      if (R > upper)
	{
	  R = upper;
	}
      else
	{
	  //setValue(R);
	  xcur=R;
	  while (k-- > 0 && pd(xcur) > z)
	    {
	      R += width;
	      if (R > upper)
		{
		  R = upper;
		  break;
		}
	      xcur = R;
	    }
	}

      T xnew;
      for(;;)
	{
	  xnew =  L + rng() * (R - L);
	  xcur=xnew;
	  T g = pd(xcur);
	  if (g >= z - std::numeric_limits<T>::epsilon())
	    {
	      break;
	    }
	  else
	    {
	      if (xnew < xold)
		{
		  L = xnew;
		}
	      else
		{
		  R = xnew;
		}
	    }
	}
      
      if (adapt)
	{
	  sumdiff += niter * std::abs(xnew - xold);
	  ++niter;
	  if ( niter > min_adapt)
	    {
	      width = 2 * sumdiff / niter / (niter - 1);  
	    }
	}
      
      return xcur;
    }

  };
}

#endif
