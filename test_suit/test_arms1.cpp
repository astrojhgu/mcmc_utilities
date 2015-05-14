#include "core/arms.hpp"
#include <algorithm>
using namespace std;
using namespace mcmc_utilities;
/* ********************************************************************* */


/* ********************************************************************* */

/* ********************************************************************* */

class normix
  :public probability_density_md<double,double>
{
private:
  normix* do_clone()const
  {
    return new normix(*this);
  }

  void do_var_range(double& x1,double& x2)const
  {
    x1=-40;
    x2=40;
  }

  double do_eval(const double& x)const
  {
    return exp(-x*x/(2*.001*.001));//+1e-199;
    double p1 = .3;
    double mean1 = 4.;
    double sd1 = .01;
    double mean2 = 10.;
    double sd2 = 2.5;

    return (p1/sd1)*std::exp(-0.5*std::pow(((x - mean1)/sd1),2.0))
      + ((1-p1)/sd2)*std::exp(-0.5*std::pow(((x - mean2)/sd2),2.0));
  }
}norm;



/* ********************************************************************* */

int main(void)
{
  FILE *f;
  int err,  npoint = 500, nsamp = 1;
  int neval,i;
  std::vector<double> xinit;
  double  xl,xr;
  norm.var_range(xl,xr);
  for(double x=xl+.1;x<xr;x+=(xr-xl)/10)
    {
      xinit.push_back(x);
    }
  for(double x=-1;x<=1;x+=.1)
    {
      xinit.push_back(x);
    }
  sort(xinit.begin(),xinit.end());
  //double xcent[10], qcent[10] = {5., 30., 70., 95.};
  
  std::vector<double> qcent;
  //qcent.push_back(5);
  //qcent.push_back(30);
  //qcent.push_back(70);
  //qcent.push_back(95);
  std::vector<double> xcent(qcent.size());
  std::vector<double> xsamp(10);
  unsigned seed = 44;
  double convex = 1.0;
  int dometrop = 1;
  double xprev = 0.0;

  /* initialise random number generator */
  srand(seed);

  /* open a file */
  f = fopen("arms.out02","w+");

  for(i=0;i<100000;i++){
    err = arms(xinit,xl,xr,norm,convex,
	       npoint,dometrop,xprev,xsamp,qcent,xcent,neval, u_random<double>);
    //err= arms_simple (300, norm, xprev,xsamp);
    //const probability_density_md<T_p,T_var>& myfunc,
    //int dometrop, T_var *xprev, T_var *xsamp);
    
    if(err>0){
      fprintf(f,"error code = %d\n",err);
      exit(1);
    }
    /* update xprev to get a Markov chain */
    xprev = xsamp[0];
    fprintf(f,"%d  %11.5f  %d\n",i,xsamp[0],neval);
  }
  return 0;
}

