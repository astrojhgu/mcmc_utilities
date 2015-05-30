#ifndef CORE_ARMS_HPP
#define CORE_ARMS_HPP

/* adaptive rejection metropolis sampling */
/*
  Downloaded from
  http://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html
  References
  Adaptive Rejection Sampling
by Wally Gilks wally.gilks@maths.leeds.ac.uk
Department of Statistics, University of Leeds, UK
Adaptive rejection sampling (ARS) is a method for efficiently sampling from any univariate probability density function which is log-concave. It is very useful in applications of Gibbs sampling, where full-conditional distributions are algebraically very messy yet often log-concave. For non-log-concave distributions, ARS can be followed by a single step of the Metropolis-Hastings algorithm; this is then adaptive rejection Metropolis sampling (ARMS).
Adaptive rejection sampling (ARS)

ARS works by constructing an envelope function of the log of the target density, which is then used in rejection sampling (see, for example, Ripley, 1987). Whenever a point is rejected by ARS, the envelope is updated to correspond more closely to the true log density, thereby reducing the chance of rejecting subsequent points. Fewer ARS rejection steps implies fewer point-evaluations of the log density. Such evaluations are typically very expensive computationally in applications of Gibbs sampling.
In the original formulation of ARS, the envelope is constructed from a set of tangents to the log-density (Gilks and Wild, 1992). In a later version the envelope is constructed from chords (secants) intersecting on the log-density (Gilks, 1992). Both methods assume that the log density is concave, which is surprisingly often true even for very messy full conditional distributions encountered in Gibbs sampling. ARS is the principle sampling methodology used in the BUGS/WinBUGS project software.

I have written some FORTRAN code for adaptive rejection sampling from log-concave densities, using the 'tangent-based' method.

Adaptive Metropolis rejection sampling (ARMS)

Introduction

Adaptive rejection Metropolis sampling (ARMS) is a method for efficiently sampling from complicated univariate densities, such as typically occur in applications of Gibbs sampling (Gilks, Best and Tan, 1995). ARMS is a generalisation of the method of adaptive rejection sampling (ARS) (Gilks, 1992) , which was itself a development of the original method proposed by Gilks and Wild (1992). The ARMS generalisation includes a Metropolis step to accomodate non-concavity in the log density.
These notes accompany C code implementing ARMS. It should be possible to call this C code from a FORTRAN program, but I haven't tried this yet.

Method

ARS works by constructing an envelope function of the log of the target density, which is then used in rejection sampling (see, for example, Ripley, 1987). Whenever a point is rejected by ARS, the envelope is updated to correspond more closely to the true log density, thereby reducing the chance of rejecting subsequent points. Fewer ARS rejection steps implies fewer point-evaluations of the log density. Such evaluations are typically very expensive computationally in applications of Gibbs sampling.

In the original formulation of ARS, the envelope is constructed from a set of tangents to the log-density (Gilks and Wild, 1992). In a later version the envelope is constructed from chords (secants) intersecting on the log-density (Gilks, 1992). Both methods assume that the log density is concave, which is surprisingly often true even for very messy full conditional distributions encountered in Gibbs sampling.

Occasionally however log-concavity does not obtain, typically in non-linear models, or with non-exponential- family distributions. This can be particularly frustrating when the target density is known to be very nearly log-concave. ARMS deals with this situation by performing a Metropolis step (Metropolis et al, 1953) on each point accepted at an ARS rejection step. In the Metropolis step, the new point is weighed against the previous point sampled (in applications of Gibbs sampling this is the previous value of the model parameter currently being updated by Gibbs). If the new point is rejected, the previous point is retained as the new point. The procedure is exact (in the sense of returning samples from the exact target density), regardless of the degree of convexity in the log density. However it is most efficient (rejecting fewer points at the Metropolis step) when the density is approximately log-concave. Indeed, when the density is truly log-concave, the Metropolis step will never reject.

There is, however, an overhead in using ARMS instead of ARS when the log density is truly log-concave. This is for two reasons. Firstly, each call to the ARMS function will require an additional point-evaluation of the log-density (to prime the system for Metropolis); and secondly, when log- concavity is assured, squeezing functions can be constructed which may save a function evaluation at each ARS-rejection- step. If sampling only one point from the density, (as is typically the case for Gibbs sampling) the loss in efficiency will be small. If many samples from the same density are required, the loss in efficiency would be severe. Therefore the ARMS function includes a parameter which allows you to indicate if the density is known to be log-concave. If indicated thus, the ARMS function will not implement the Metropolis step, and will exit if non-concavity is detected.

The ARMS code allows you to choose initial construction points for the envelope; (I will continue to call the ARMS rejection function an 'envelope', although for non-log-concave distributions this function is not in general an envelope any more). You can set these initial construction points anywhere you like; however, in applications of Gibbs sampling, when the distribution to be sampled is not log-concave, it is important that initial construction points are set independently of the current value of the parameter being updated. (CAUTION: the advice given in Gilks, Best and Tan, 1995, to set initial construction points based on the previous Gibbs iterate, is only valid if the target density is log-concave; a corrigendum has been submitted to Applied Statistics). In many applications, for each variable to be updated by ARMS, a fixed set of 4 or 5 initial construction points will suffice. When this would result in a large number of iterations of ARS, some exploration of the density to be sampled may assist in choosing an efficient set of initial construction points.

Both ARS and ARMS can be used to generate many points from the target density (although in applications of Gibbs sampling only one point need be sampled from each target density). As the envelope function used by both methods progressively adapts to the shape of the target density, sampling becomes progressively more efficient as more points are sampled from the target density. If the target density is log-concave, the sampled points will be independent from the exact target density. If there is non-log- concavity, the samples will only be 'from' the target density in the limit, i.e. only after many samples have been drawn, and even then they will be dependent samples. Note that this is not a concern in applications of Gibbs sampling: only one point need be sampled as the stationary distribution (of the induced Markov chain) is preserved by ARMS.


  Gilks, W. R. (1992) Derivative-free adaptive rejection sampling for Gibbs sampling. Bayesian Statistics 4, (eds. Bernardo, J., Berger, J., Dawid, A. P., and Smith, A. F. M.) Oxford University Press.
  Gilks, W. R., Best, N. G. and Tan, K. K. C. (1995) Adaptive rejection Metropolis sampling. Applied Statistics, 44, 455-472.
  
  Gilks, W. R. and Wild, P. (1992) Adaptive rejection sampling for Gibbs sampling. Applied Statistics 41, pp 337-348.
  
  Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H. and Teller, E. (1953) Equations of state calculations by fast computing machines. J. Chem. Phys., 21, 1087-1092.
  
  Ripley, B. (1987) Stochastic Simulation. New York, Wiley. 
*/
/* *********************************************************************** */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <vector>
#include "distribution.hpp"
#include "mcmc_exception.hpp"
#include "find_peak.hpp"
/* *********************************************************************** */

namespace mcmc_utilities
{
  template <typename T_p,typename T_var>
  struct point {    /* a point in the x,y plane */
    T_var x;
    T_p y;             /* x and y coordinates */
    T_p ey;              /* exp(y-ymax+YCEIL) */
    T_p cum;             /* integral up to x of rejection envelope */
    int f;                  /* is y an evaluated point of log-density */
    point *pl,*pr;   /* envelope points to left and right of x */
  };
  
  /* *********************************************************************** */
  
  template <typename T_p,typename T_var>
  struct envelope {  /* attributes of the entire rejection envelope */
    int cpoint;              /* number of POINTs in current envelope */
    int npoint;              /* max number of POINTs allowed in envelope */
    int *neval;              /* number of function evaluations performed */
    T_p ymax;             /* the maximum y-value in the current envelope */
    point<T_p,T_var> *p;                /* start of storage of envelope POINTs */
    const T_p *convex;          /* adjustment for convexity */

    envelope()
      :p(0)
    {}

    ~envelope()
    {
      if(p)
	{
	  delete[] p;
	}
    }
  };
  
  /* *********************************************************************** */
  
  /* *********************************************************************** */
  
  template <typename T_p,typename T_var>
  struct metropolis { /* for metropolis step */
    bool on;            /* whether metropolis is to be used */
    T_var xprev;      /* previous Markov chain iterate */
    T_p yprev;      /* current log density at xprev */
  };
  
  /* *********************************************************************** */
  
  //#define RAND_MAX 2147483647      /* For Sun4 : remove this for some systems */
  template <typename T>
  struct var_constants
  {
    static const T XEPS()
    {
      return .00001;
    }
  };/* critical relative x-value difference */

  template <typename T>
  struct prob_constants
  {
    static const T YEPS()
    {
      return .01;
    }                /* critical y-value difference */
    static const T EYEPS()
    {
      return .001;
    }              /* critical relative exp(y) difference */
    static const T YCEIL()
    {
      //return std::log(std::numeric_limits<T>::max()/10000.);
      return 100;
    }                /* maximum y avoiding overflow in exp(y) */
  };
  /* *********************************************************************** */
  
  /* declarations for functions defined in this file */
  template <typename T_p,typename T_var, typename T_urand>
  int arms_simple (int ninit,const probability_density_1d<T_p,T_var>& myfunc,
		   const T_var& xprev, std::vector<T_var>& xsamp, bool dometrop, const T_urand& urand,bool use_peak_finder=true);

  
  template <typename T_p,typename T_var, typename T_urand>
  int arms (const std::vector<T_var>& xinit, const T_var& xl, const T_var& xr, 
	    const probability_density_1d<T_p,T_var>& myfunc,
	    const T_p& convex, int npoint, bool dometrop, const T_var& xprev, std::vector<T_var>& xsamp,
	    int& neval,const T_urand& urand);

  template <typename T_p,typename T_var>
  int initial (const std::vector<T_var>& xinit, T_var xl, T_var xr, int npoint,
	       const probability_density_1d<T_p,T_var>& myfunc, envelope<T_p,T_var>& env, const T_p& convex, int& neval,
	       metropolis<T_p,T_var>& metrop);
  template <typename T_p,typename T_var,typename T_urand>
  void sample(envelope<T_p,T_var>& env, point<T_p,T_var> *p,const T_urand& urand);
  
  template <typename T_p,typename T_var>
  void invert(T_p prob, envelope<T_p,T_var>& env, point<T_p,T_var> *p);

  template <typename T_p,typename T_var,typename T_urand>
  int test(envelope<T_p,T_var>& env, point<T_p,T_var> *p, const probability_density_1d<T_p,T_var>& myfunc, metropolis<T_p,T_var>& metrop, const T_urand& urand);
  
  template <typename T_p,typename T_var>
  int update(envelope<T_p,T_var>& env, point<T_p,T_var> *p, const probability_density_1d<T_p,T_var>& myfunc, metropolis<T_p,T_var>& metrop);

  template <typename T_p,typename T_var>
  void cumulate(envelope<T_p,T_var>& env);

  template <typename T_p,typename T_var>
  int meet (point<T_p,T_var> *q, envelope<T_p,T_var>& env, metropolis<T_p,T_var>& metrop);

  template <typename T_p,typename T_var>
  T_p area(point<T_p,T_var> *q);


  template <typename T_p>
  T_p expshift(T_p y, T_p y0);

  template <typename T_p>
  T_p logshift(T_p y, T_p y0);


  template <typename T_p,typename T_var>
  T_p perfunc(const probability_density_1d<T_p,T_var>& myfunc, envelope<T_p,T_var>& env, T_var x);

  template <typename T_var>
  class u_random;

  /* *********************************************************************** */

  template <typename T_p,typename T_var, typename T_urand>
  int arms_simple (int ninit,const probability_density_1d<T_p,T_var>& myfunc,
		   const T_var& xprev, std::vector<T_var>& xsamp, bool dometrop, const T_urand& urand,bool use_peak_finder)
    
  /* adaptive rejection metropolis sampling - simplified argument list */
  /* ninit        : number of starting values to be used */
  /* *xl          : left bound */
  /* *xr          : right bound */
  /* *myfunc      : function to evaluate log-density */
  /* *mydata      : data required by *myfunc */
  /* dometrop     : whether metropolis step is required */
  /* *xprev       : current value from markov chain */
  /* *xsamp       : to store sampled value */
    
  {
    std::vector<T_var> xinit(ninit);
    T_var xl,xr;
    
    myfunc.var_range(xl,xr);
    if(xl>=xr)
      {
	throw range_not_ordered();
      }
    if(use_peak_finder)
      {
	for(int i=0;i<ninit-1;++i)
	  {
	    xinit[i]=xl+(xr-xl)/(ninit+1)*(i+1);
	    //xinit[i]=xl_shifted+(xr_shifted-xl_shifted)/(ninit+2)*(i+1);
	  }
	xinit[ninit-1]=find_peak(myfunc);
	std::sort(xinit.begin(),xinit.end());
      }
    else
      {
	for(int i=0;i<ninit;++i)
	  {
	    xinit[i]=xl+(xr-xl)/(ninit+1)*(i+1);
	    //xinit[i]=xl_shifted+(xr_shifted-xl_shifted)/(ninit+2)*(i+1);
	  }
      }
    int npoint=std::max(200,2*ninit + 2);
    T_var convex=1.;
    //int dometrop=1;
    int err=0;
    int neval;
    err = arms(xinit,xl,xr,myfunc,convex,
	       npoint,dometrop,xprev,xsamp,neval, urand);
    
    return err;
  }

  /* *********************************************************************** */
  
  template <typename T_p,typename T_var, typename T_urand>
  int arms (const std::vector<T_var>& xinit, const T_var& xl, const T_var& xr, 
	    const probability_density_1d<T_p,T_var>& myfunc,
	    const T_p& convex, int npoint, bool dometrop, const T_var& xprev, std::vector<T_var>& xsamp,
	    int& neval,const T_urand& urand)
    
  /* to perform derivative-free adaptive rejection sampling with metropolis step */
  /* *xinit       : starting values for x in ascending order */
  /* ninit        : number of starting values supplied */
  /* *xl          : left bound */
  /* *xr          : right bound */
  /* *myfunc      : function to evaluate log-density */
  /* *mydata      : data required by *myfunc */
  /* *convex      : adjustment for convexity */
  /* npoint       : maximum number of envelope points */
  /* dometrop     : whether metropolis step is required */
  /* *xprev       : previous value from markov chain */
  /* *xsamp       : to store sampled values */
  /* nsamp        : number of sampled values to be obtained */
  /* *qcent       : percentages for envelope centiles */
  /* *xcent       : to store requested centiles */
  /* ncent        : number of centiles requested */
  /* *neval       : on exit, the number of function evaluations performed */
    
  {
    //int ninit=xinit.size();
    int nsamp=xsamp.size();
    //envelope<T_p,T_var> *env;      /* rejection envelope */
    envelope<T_p,T_var> env;
    point<T_p,T_var> pwork;        /* a working point, not yet incorporated in envelope */
    int msamp=0;        /* the number of x-values currently sampled */
    //funbag<T_p,T_var> lpdf;        /* to hold density function and its data */
    metropolis<T_p,T_var> metrop; /* to hold bits for metropolis step */
    int i,err;
    
    
    /* incorporate density function and its data into funbag lpdf */
    //lpdf.mydata = mydata;
    //lpdf.myfunc = myfunc;
    
    /* set up space required for envelope */
    //env = (envelope *)malloc(sizeof(envelope));
    //env=new envelope<T_p,T_var>;
    //if(env == NULL){
    //  /* insufficient space */
    //  return 1006;
    //}
    
    /* start setting up metropolis struct */
    //metrop = (metropolis *)malloc(sizeof(metropolis));
    //metrop=new metropolis<T_p,T_var>;
    //if(metrop == NULL){
      /* insufficient space */
    //  return 1006;
    //}
    metrop.on = dometrop;
    
    /* set up initial envelope */
    err = initial(xinit,xl,xr,npoint,myfunc,env,convex,
		  neval,metrop);
    //std::ofstream ofs("debug.txt");
    //ofs<<env;
    //exit(0);
    if(err)return err;
    
    /* finish setting up metropolis struct (can only do this after */
    /* setting up env) */
    if(metrop.on){
      if((xprev < xl) || (xprev > xr)){
	/* previous markov chain iterate out of range */
	//return 1007;
	throw previous_iter_out_of_range();
      }
      metrop.xprev = xprev;
      metrop.yprev = perfunc(myfunc,env,xprev);
    }
    
    /* now do adaptive rejection */
    do {
      /* sample a new point */
      sample (env,&pwork, urand);
      
      /* perform rejection (and perhaps metropolis) tests */
      i = test(env,&pwork,myfunc,metrop,urand);
      if(i == 1){
	/* point accepted */
	xsamp[msamp++] = pwork.x;
      } else if (i != 0) {
	/* envelope error - violation without metropolis */
	//return 2000;
	throw envelope_error();
      }  
    } while (msamp < nsamp);
    
    /* free space */
    //free(env->p);
    //free(env);
    //free(metrop);
    //delete[] env.p;
    //delete env;
    //delete metrop;
    
    return 0;
  }
  
  /* *********************************************************************** */
  
  template <typename T_p,typename T_var>
  int initial (const std::vector<T_var>& xinit, T_var xl, T_var xr, int npoint,
	       const probability_density_1d<T_p,T_var>& myfunc, envelope<T_p,T_var>& env, const T_p& convex, int& neval,
	       metropolis<T_p,T_var>& metrop)
    
  /* to set up initial envelope */
  /* xinit        : initial x-values */
  /* ninit        : number of initial x-values */
  /* xl,xr        : lower and upper x-bounds */
  /* npoint       : maximum number of POINTs allowed in envelope */
  /* *lpdf        : to evaluate log density */
  /* *env         : rejection envelope attributes */
  /* *convex      : adjustment for convexity */
  /* *neval       : current number of function evaluations */
  /* *metrop      : for metropolis step */
    
  {
    int ninit=xinit.size();
    int i,j,k,mpoint;
    point<T_p,T_var> *q;
    
    if(ninit<3){
      /* too few initial points */
      //return 1001;
      throw too_few_init_points();
    }
    
    mpoint = 2*ninit + 1;
    if(npoint < mpoint){
      /* too many initial points */
      //return 1002;
      throw too_many_init_points();
    }
    
    if((xinit[0] <= xl) || (xinit[ninit-1] >= xr)){
      /* initial points do not satisfy bounds */
      //return 1003;
      throw init_point_out_of_range();
    }
    
    for(i=1; i<ninit; i++){
      if(xinit[i] <= xinit[i-1]){
	/* data not ordered */
	//return 1004;
	throw data_not_ordered();
      }
    }
    
    if(convex < 0.0){
      /* negative convexity parameter */
      ///return 1008;
      
    }
    
    /* copy convexity address to env */
    env.convex = &convex;
    
    /* copy address for current number of function evaluations */
    env.neval = &neval;
    /* initialise current number of function evaluations */
    *(env.neval) = 0;
    
    /* set up space for envelope POINTs */
    env.npoint = npoint;
    //env->p = (point *)malloc(npoint*sizeof(point));
    env.p = new point<T_p,T_var>[npoint];
    if(env.p == NULL){
      /* insufficient space */
      //return 1006;
      throw insufficient_space();
    }
    
    /* set up envelope POINTs */
    q = env.p;
    /* left bound */
    q->x = xl;
    q->f = 0;
    q->pl = NULL;
    q->pr = q+1;
    for(j=1, k=0; j<mpoint-1; j++){
      q++;
      if(j%2){
	/* point on log density */
	q->x = xinit[k++];
	q->y = perfunc(myfunc,env,q->x);
	q->f = 1;
      } else {
	/* intersection point */
	q->f = 0;
      }
      q->pl = q-1;
      q->pr = q+1;
    }
    /* right bound */
    q++;
    q->x = xr;
    q->f = 0;
    q->pl = q-1;
    q->pr = NULL;
    
    /* calculate intersection points */
    //a
    q = env.p;
    for (j=0; j<mpoint; j=j+2, q=q+2){
      //assert(!isinf(q->x));
      //cout<<q->x<<endl;
      if(meet(q,env,metrop)){
	/* envelope violation without metropolis */
	//return 2000;
	throw envelope_violation_without_metropolis();
      }
    }
    
    /* exponentiate and integrate envelope */
    cumulate(env);
    
    /* note number of points currently in envelope */
    env.cpoint = mpoint;
    
    return 0;
  }

  /* *********************************************************************** */

  template <typename T_p,typename T_var,typename T_urand>
  void sample(envelope<T_p,T_var>& env, point<T_p,T_var> *p,const T_urand& urand)

  /* To sample from piecewise exponential envelope */
  /* *env    : envelope attributes */
  /* *p      : a working point to hold the sampled value */

  {
    T_p prob;
  
    /* sample a uniform */
    //prob = u_random<T_var>();
    prob = urand();
    /* get x-value correponding to a cumulative probability prob */
    invert(prob,env,p);

    return;
  }

  /* *********************************************************************** */

  template <typename T_p,typename T_var>
  void invert(T_p prob, envelope<T_p,T_var>& env, point<T_p,T_var> *p)

  /* to obtain a point corresponding to a qiven cumulative probability */
  /* prob    : cumulative probability under envelope */
  /* *env    : envelope attributes */
  /* *p      : a working point to hold the sampled value */

  {
    T_var xl(0),xr(0);
    T_p u(0),yl(0),yr(0),eyl(0),eyr(0),prop(0);
    point<T_p,T_var> *q;

    /* find rightmost point in envelope */
    q = env.p;
    while(q->pr != NULL)q = q->pr;
    
    /* find exponential piece containing point implied by prob */
    u = prob * q->cum;
    while(q->pl->cum > u)
      {
	q = q->pl;
      }
    /* piece found: set left and right points of p, etc. */
    p->pl = q->pl;
    p->pr = q;
    p->f = 0;
    p->cum = u;

    /* calculate proportion of way through integral within this piece */
    prop = (u - q->pl->cum) / (q->cum - q->pl->cum);

    /* get the required x-value */
    if (q->pl->x == q->x){
      /* interval is of zero length */
      p->x = q->x;
      p->y = q->y;
      p->ey = q->ey;
    } else {
      xl = q->pl->x;
      xr = q->x;
      yl = q->pl->y;
      yr = q->y;
      eyl = q->pl->ey;
      eyr = q->ey;
      if(std::abs(yr - yl) < prob_constants<T_p>::YEPS()){
	/* linear approximation was used in integration in function cumulate */
	if(std::abs(eyr - eyl) > prob_constants<T_p>::EYEPS()*std::abs(eyr + eyl)){
	  p->x = xl + ((xr - xl)/(eyr - eyl))
	    * (-eyl + sqrt((1. - prop)*eyl*eyl + prop*eyr*eyr));
	} else {
	  p->x = xl + (xr - xl)*prop;
	}
	p->ey = ((p->x - xl)/(xr - xl)) * (eyr - eyl) + eyl;
	p->y = logshift(p->ey, env.ymax);
      } else {
	/* piece was integrated exactly in function cumulate */
	p->x = xl + ((xr - xl)/(yr - yl))
	  * (-yl + logshift(((1.-prop)*eyl + prop*eyr), env.ymax));
	p->y = ((p->x - xl)/(xr - xl)) * (yr - yl) + yl;
	p->ey = expshift(p->y, env.ymax);
      }
    }

    p->x=std::min(p->x,xr);
    p->x=std::max(p->x,xl);
    /* guard against imprecision yielding point outside interval */
    if ((p->x < xl) || (p->x > xr))
      {
	throw arms_exception(3);
      };

    return;
  }

  /* *********************************************************************** */

  template <typename T_p,typename T_var,typename T_urand>
  int test(envelope<T_p,T_var>& env, point<T_p,T_var> *p, const probability_density_1d<T_p,T_var>& myfunc, metropolis<T_p,T_var>& metrop, const T_urand& urand)

  /* to perform rejection, squeezing, and metropolis tests */
  /* *env          : envelope attributes */
  /* *p            : point to be tested */
  /* *lpdf         : to evaluate log-density */
  /* *metrop       : data required for metropolis step */

  {
    T_p u,y,ysqueez,ynew,yold,znew,zold,w;
    point<T_p,T_var> *ql,*qr;
  
    /* for rejection test */
    //u = u_random<T_var>() * p->ey;
    u = urand() * p->ey;
    y = logshift(u,env.ymax);

    if(!(metrop.on) && (p->pl->pl != NULL) && (p->pr->pr != NULL)){
      /* perform squeezing test */
      if(p->pl->f){
	ql = p->pl;
      } else {
	ql = p->pl->pl;
      }
      if(p->pr->f){
	qr = p->pr;
      } else {
	qr = p->pr->pr;
      }
      ysqueez = (qr->y * (p->x - ql->x) + ql->y * (qr->x - p->x))
	/(qr->x - ql->x);
      if(y <= ysqueez){
	/* accept point at squeezing step */
	return 1;
      }
    }

    /* evaluate log density at point to be tested */
    ynew = perfunc(myfunc,env,p->x);
  
    /* perform rejection test */
    if(!(metrop.on) || ((metrop.on) && (y >= ynew))){
      /* update envelope */
      p->y = ynew;
      p->ey = expshift(p->y,env.ymax);
      p->f = 1;
      if(update(env,p,myfunc,metrop)){
	/* envelope violation without metropolis */
	//return -1;
	throw envelope_violation_without_metropolis();
      }
      /* perform rejection test */
      if(y >= ynew){
	/* reject point at rejection step */
	return 0;
      } else {
	/* accept point at rejection step */
	return 1;
      }
    }

    /* continue with metropolis step */
    yold = metrop.yprev;
    /* find envelope piece containing metrop->xprev */
    ql = env.p;
    while(ql->pl != NULL)ql = ql->pl;
    while(ql->pr->x < metrop.xprev)ql = ql->pr;
    qr = ql->pr;
    /* calculate height of envelope at metrop->xprev */
    w = (metrop.xprev - ql->x)/(qr->x - ql->x);
    zold = ql->y + w*(qr->y - ql->y);
    znew = p->y;
    if(yold < zold)zold = yold;
    if(ynew < znew)znew = ynew;
    w = ynew-znew-yold+zold;
    if(w > 0.0)w = 0.0;

    if(w > -prob_constants<T_p>::YCEIL()){
      w = std::exp(w);
    } else {
      w = 0.0;
    }
    //u = u_random<T_var>();
    u = urand();
    if(u > w){
      /* metropolis says dont move, so replace current point with previous */
      /* markov chain iterate */
      p->x = metrop.xprev;
      p->y = metrop.yprev;
      p->ey = expshift(p->y,env.ymax);
      p->f = 1;
      p->pl = ql;
      p->pr = qr;
    } else {
      /* trial point accepted by metropolis, so update previous markov */
      /* chain iterate */
      metrop.xprev = p->x;
      metrop.yprev = ynew;
    }
    return 1;
  }

  /* *********************************************************************** */

  template <typename T_p,typename T_var>
  int update(envelope<T_p,T_var>& env, point<T_p,T_var> *p, const probability_density_1d<T_p,T_var>& myfunc, metropolis<T_p,T_var>& metrop)

  /* to update envelope to incorporate new point on log density*/
  /* *env          : envelope attributes */
  /* *p            : point to be incorporated */
  /* *lpdf         : to evaluate log-density */
  /* *metrop       : for metropolis step */

  {
    point<T_p,T_var> *m,*ql,*qr,*q;

    if(!(p->f) || (env.cpoint > env.npoint - 2)){
      /* y-value has not been evaluated or no room for further points */
      /* ignore this point */
      return 0;
    }

    /* copy working point p to a new point q */
    q = env.p + env.cpoint++;
    q->x = p->x;
    q->y = p->y;
    q->f = 1;

    /* allocate an unused point for a new intersection */
    m = env.p + env.cpoint++;
    m->f = 0;
    if((p->pl->f) && !(p->pr->f)){
      /* left end of piece is on log density; right end is not */
      /* set up new intersection in interval between p->pl and p */
      m->pl = p->pl;
      m->pr = q;
      q->pl = m;
      q->pr = p->pr;
      m->pl->pr = m;
      q->pr->pl = q;
    } else if (!(p->pl->f) && (p->pr->f)){
      /* left end of interval is not on log density; right end is */
      /* set up new intersection in interval between p and p->pr */
      m->pr = p->pr;
      m->pl = q;
      q->pr = m;
      q->pl = p->pl;
      m->pr->pl = m;
      q->pl->pr = q;
    } else {
      /* this should be impossible */
      throw arms_exception(10);
    }

    /* now adjust position of q within interval if too close to an endpoint */
    if(q->pl->pl != NULL){
      ql = q->pl->pl;
    } else {
      ql = q->pl;
    }
    if(q->pr->pr != NULL){
      qr = q->pr->pr;
    } else {
      qr = q->pr;
    }
    if (q->x < (1. - var_constants<T_var>::XEPS()) * ql->x + var_constants<T_var>::XEPS() * qr->x){
      /* q too close to left end of interval */
      q->x = (1. - var_constants<T_var>::XEPS()) * ql->x + var_constants<T_var>::XEPS() * qr->x;
      q->y = perfunc(myfunc,env,q->x);
    } else if (q->x > var_constants<T_var>::XEPS() * ql->x + (1. - var_constants<T_var>::XEPS()) * qr->x){
      /* q too close to right end of interval */
      q->x = var_constants<T_var>::XEPS() * ql->x + (1. - var_constants<T_var>::XEPS()) * qr->x;
      q->y = perfunc(myfunc,env,q->x);
    }

    /* revise intersection points */
    if(meet(q->pl,env,metrop)){
      /* envelope violation without metropolis */
      //return 1;
      throw envelope_violation_without_metropolis();
    }
    if(meet(q->pr,env,metrop)){
      /* envelope violation without metropolis */
      //return 1;
      throw envelope_violation_without_metropolis();
    }
    if(q->pl->pl != NULL){
      if(meet(q->pl->pl->pl,env,metrop)){
	/* envelope violation without metropolis */
	//return 1;
	throw envelope_violation_without_metropolis();
      }
    }
    if(q->pr->pr != NULL){
      if(meet(q->pr->pr->pr,env,metrop)){
	/* envelope violation without metropolis */
	//return 1;
	throw envelope_violation_without_metropolis();
      }
    }

    /* exponentiate and integrate new envelope */
    cumulate(env);

    return 0;
  }

  /* *********************************************************************** */

  template <typename T_p,typename T_var>
  void cumulate(envelope<T_p,T_var>& env)

  /* to exponentiate and integrate envelope */
  /* *env     : envelope attributes */

  {
    point<T_p,T_var> *q,*qlmost;

    qlmost = env.p;
    /* find left end of envelope */
    while(qlmost->pl != NULL)qlmost = qlmost->pl;

    /* find maximum y-value: search envelope */
    env.ymax = qlmost->y;
    for(q = qlmost->pr; q != NULL; q = q->pr){
      if(q->y > env.ymax)env.ymax = q->y;
    }

    /* exponentiate envelope */
    for(q = qlmost; q != NULL; q = q->pr){
      q->ey = expshift(q->y,env.ymax);
    }

    /* integrate exponentiated envelope */
    qlmost->cum = 0.;
    for(q = qlmost->pr; q != NULL; q = q->pr){
      q->cum = q->pl->cum + area(q);
    }

    return;
  }

  /* *********************************************************************** */

  template <typename T_p,typename T_var>
  int meet (point<T_p,T_var> *q, envelope<T_p,T_var>& env, metropolis<T_p,T_var>& metrop)
  /* To find where two chords intersect */
  /* q         : to store point of intersection */
  /* *env      : envelope attributes */
  /* *metrop   : for metropolis step */

  {
    T_p gl(0),gr(0),grl(0),dl(0),dr(0);
    int il,ir,irl;

    if(q->f){
      /* this is not an intersection point */
      throw arms_exception(30);
    }

    /* calculate coordinates of point of intersection */
    if ((q->pl != NULL) && (q->pl->pl->pl != NULL)){
      /* chord gradient can be calculated at left end of interval */
      gl = (q->pl->y - q->pl->pl->pl->y)/(q->pl->x - q->pl->pl->pl->x);
      //cout<<"d:"<<q->pl->y<<" "<<q->pl->pl->pl->y<<endl;
      il = 1;
    } else {
      /* no chord gradient on left */
      il = 0;
    }
    if ((q->pr != NULL) && (q->pr->pr->pr != NULL)){
      /* chord gradient can be calculated at right end of interval */
      gr = (q->pr->y - q->pr->pr->pr->y)/(q->pr->x - q->pr->pr->pr->x);
      ir = 1;
    } else {
      /* no chord gradient on right */
      ir = 0;
    }
    if ((q->pl != NULL) && (q->pr != NULL)){
      /* chord gradient can be calculated across interval */
      grl = (q->pr->y - q->pl->y)/(q->pr->x - q->pl->x);
      irl = 1;
    } else {
      irl = 0;
    }

    if(irl && il && (gl<grl)){
      /* convexity on left exceeds current threshold */
      if(!(metrop.on)){
	/* envelope violation without metropolis */
	//return 1;
	throw envelope_violation_without_metropolis();
      }
      /* adjust left gradient */
      gl = gl + (1.0 + *(env.convex)) * (grl - gl);
    }

    if(irl && ir && (gr>grl)){
      /* convexity on right exceeds current threshold */
      if(!(metrop.on)){
	/* envelope violation without metropolis */
	//return 1;
	throw envelope_violation_without_metropolis();
      }
      /* adjust right gradient */
      gr = gr + (1.0 + *(env.convex)) * (grl - gr);
    }

    if(il && irl){
      dr = (gl - grl) * (q->pr->x - q->pl->x);
      if(dr < prob_constants<T_p>::YEPS()){
	/* adjust dr to avoid numerical problems */
	dr = prob_constants<T_p>::YEPS();
      }
    }

    if(ir && irl){
      dl = (grl - gr) * (q->pr->x - q->pl->x);
      if(dl < prob_constants<T_p>::YEPS()){
	/* adjust dl to avoid numerical problems */
	dl = prob_constants<T_p>::YEPS();
      }
    }

    if(il && ir && irl){
      /* gradients on both sides */
      q->x = (dl * q->pr->x + dr * q->pl->x)/(dl + dr);
      //cout<<q->x<<endl;
      //cout<<"a:"<<q->pr->x<<" "<<q->pl->x<<" "<<q->x<<endl;
      //cout<<"b:"<<dl<<" "<<dr<<endl;
      //cout<<"c:"<<gl<<" "<<gr<<endl;
      q->y = (dl * q->pr->y + dr * q->pl->y + dl * dr)/(dl + dr);
    } else if (il && irl){
      /* gradient only on left side, but not right hand bound */
      q->x = q->pr->x;
      q->y = q->pr->y + dr;
    } else if (ir && irl){
      /* gradient only on right side, but not left hand bound */
      q->x = q->pl->x;
      q->y = q->pl->y + dl;
    } else if (il){
      /* right hand bound */
      q->y = q->pl->y + gl * (q->x - q->pl->x);
    } else if (ir){
      /* left hand bound */
      q->y = q->pr->y - gr * (q->pr->x - q->x);
    } else {
      /* gradient on neither side - should be impossible */
      throw arms_exception(31);
    }
    //cerr<<q->x<<endl;
    if(((q->pl != NULL) && (q->x < q->pl->x)) ||
       ((q->pr != NULL) && (q->x > q->pr->x))){
      /* intersection point outside interval (through imprecision) */
      throw arms_exception(32);
      //exit(32);
    }
    /* successful exit : intersection has been calculated */
    return 0;
  }

  /* *********************************************************************** */

  template <typename T_p,typename T_var>
  T_p area(point<T_p,T_var> *q)

  /* To integrate piece of exponentiated envelope to left of point q */

  {
    T_p a;

    if(q->pl == NULL){
      /* this is leftmost point in envelope */
      //exit(1);
      throw arms_exception(2);
    } else if(q->pl->x == q->x){
      /* interval is zero length */
      a = 0.;
    } else if (std::abs(q->y - q->pl->y) < prob_constants<T_p>::YEPS()){
      /* integrate straight line piece */
      a = 0.5*(q->ey + q->pl->ey)*(q->x - q->pl->x);
    } else {
      /* integrate exponential piece */
      a = ((q->ey - q->pl->ey)/(q->y - q->pl->y))*(q->x - q->pl->x);
    }
    return a;
  }

  /* *********************************************************************** */

  template <typename T_p>
  T_p expshift(T_p y, T_p y0)

  /* to exponentiate shifted y without underflow */
  {
    if(y - y0 > -2.0 * prob_constants<T_p>::YCEIL()){
      return std::exp(y - y0 + prob_constants<T_p>::YCEIL());
    } else {
      return 0.0;
    }
  }

  /* *********************************************************************** */

  template <typename T_p>
  T_p logshift(T_p y, T_p y0)

  /* inverse of function expshift */
  {
    return (std::log(y) + y0 - prob_constants<T_p>::YCEIL());
  }

  /* *********************************************************************** */

  template <typename T_p,typename T_var>
  T_p perfunc(const probability_density_1d<T_p,T_var>& myfunc, envelope<T_p,T_var>& env, T_var x)

  /* to evaluate log density and increment count of evaluations */

  /* *lpdf   : structure containing pointers to log-density function and data */
  /* *env    : envelope attributes */
  /* x       : point at which to evaluate log density */

  {
    T_p y;

    /* evaluate density function */
    //y = (lpdf->myfunc)(x,lpdf->mydata);
    y=myfunc.eval_log(x);
    if(std::isinf(y)||std::isnan(y))
      {
	//std::cerr<<x<<std::endl;
	throw invalid_prob_value();
      }
    /* increment count of function evaluations */
    (*(env.neval))++;

    return y;
  }

  /* *********************************************************************** */

  template <typename T_p,typename T_var>
  //void display(std::ostream& os, envelope<T_p,T_var> *env)
  std::ostream& operator<<(std::ostream& os,envelope<T_p,T_var>& env)

  /* to display envelope - for debugging only */
  {
    point<T_p,T_var> *q;
  
    /* print envelope attributes */
    os<<"========================================================\n"
      <<"envelope attributes:\n"
      <<"points in use ="<<env.cpoint<<" points available = "<<env.npoint
      <<"function evaluations = "<<*(env.neval)<<std::endl
      <<"ymax = "<<env.ymax<<", p = "<<env.p<<std::endl
      <<"convexity adjustment = "<<*(env.convex)<<std::endl
      <<"--------------------------------------------------------"<<std::endl;

    /* find leftmost point */
    q = env.p;
    while(q->pl != NULL)q = q->pl;
  
    /* now print each point from left to right */
    for(q = env.p; q != NULL; q = q->pr){
      os<<"point at "<<q<<", left at "<<q->pl<<", right at "<<q->pr<<std::endl;
      os<<"x = "<<q->x<<", y = "<<q->y<<", ey = "<<q->ey<<", cum = "<<q->cum<<", f = "<<q->f<<std::endl;
    }
    os<<"========================================================\n";
  
    return os;
  }
  
  /* *********************************************************************** */
  template <typename T_var>
  class u_random
  {
  public:
    u_random()
    {
      srand(time(0));
    }
  public:
    T_var operator()()const
    {
      return ((T_var)std::rand() + 0.5)/((T_var)RAND_MAX + 1.0);
    }
  };
}
/* *********************************************************************** */


#endif

