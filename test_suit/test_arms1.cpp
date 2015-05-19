#include "core/arms.hpp"
#include <algorithm>
using namespace std;
using namespace mcmc_utilities;
/* ********************************************************************* */


/* ********************************************************************* */

/* ********************************************************************* */

class normix
  :public probability_density_1d<double,double>
{
private:
  void do_var_range(double& x1,double& x2)const
  {
    x1=-40e10;
    x2=40e10;
  }

  double do_eval_log(const double& x)const
  {
    double xmean=1e10;
    return (-(x-xmean)*(x-xmean)/(2*.001*.001*1e5*1e5));//+1e-199;
  }
}norm;



/* ********************************************************************* */

int main(void)
{
  FILE *f;
  int err, nsamp = 1;
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
  std::vector<double> xsamp(10);
  unsigned seed = 44;
  double convex = 1.0;
  int dometrop = 1;
  double xprev = 0.0;

  /* initialise random number generator */
  srand(seed);

  /* open a file */
  f = fopen("arms.out02","w+");
  int npoint=100000;
  std::vector<double> data(npoint);
  for(i=0;i<100000;i++){
    err = arms(xinit,xl,xr,norm,convex,
	       npoint,dometrop,xprev,xsamp,neval, u_random<double>);
    //err= arms_simple (300, norm, xprev,xsamp);
    //const probability_density_md<T_p,T_var>& myfunc,
    //int dometrop, T_var *xprev, T_var *xsamp);
    
    if(err>0){
      fprintf(f,"error code = %d\n",err);
      exit(1);
    }
    /* update xprev to get a Markov chain */
    xprev = xsamp.back();
    //fprintf(f,"%d  %11.5f  %d\n",i,xsamp[0],neval);
    //cout<<xprev<<endl;
    data[i]=xprev;
  }
  std::sort(data.begin(),data.end());

  for(int i=0;i<npoint;++i)
    {
      cout<<data[i]-1e10<<" "<<i<<endl;
    }
  return 0;
}

