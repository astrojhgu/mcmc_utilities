#include <iostream>
#include <vector>
#include <fstream>
#include <fio.h>
#include "core/gibbs_sampler.hpp"
#include "math/distributions.hpp"
using namespace mcmc_utilities;
using namespace std;

class variable
  :public std::vector<double>
{
public:
  double& r0;
  double& x0;
  double& y0;
  double& epsilon;
  double& theta;
  double& ampl;
  double& beta;
  double& bkg;
  
  variable()
    :vector(8),
     r0(at(0)),
     x0(at(1)),
     y0(at(2)),
     epsilon(at(3)),
     theta(at(4)),
     ampl(at(5)),
     beta(at(6)),
     bkg(at(7))
  {
    r0=27;
    x0=400;
    y0=421;
    epsilon=0;
    theta=0;
    ampl=7;
    beta=.6;
    bkg=0;

  }

private:
  variable(const variable&);
  variable& operator=(const variable&);
};

class prob
  :public probability_density_md<double,variable>
{
private:
  blitz::Array<double,2> img;
  blitz::Array<double,2> expmap;
public:
  prob(const char* img_name,const char* exp_name)
  {
    cfitsfile ff;
    ff.open(img_name);
    ff>>img;
    ff.close();
    ff.open(exp_name);
    ff>>expmap;
    ff.close();

  }
  
  double do_eval_log(const variable& param)const
  {
    static int cnt=0;
    ++cnt;
    //cout<<cnt<<endl;
    double logp=0;
#pragma omp parallel for reduction(+:logp)
    for(int i=0;i<img.extent(0);++i)
      {
	for(int j=0;j<img.extent(1);++j)
	  {
	    if(expmap(i,j)<=0)
	      {
		continue;
	      }
	       
	    double x=j;
	    double y=i;
	    double x_new=(x-param.x0)*cos(param.theta)+(y-param.y0)*sin(param.theta);
	    double y_new=(y-param.y0)*cos(param.theta)-(x-param.x0)*sin(param.theta);
	    double r_sq=x_new*x_new/exp(param.epsilon)+
	      y_new*y_new*exp(param.epsilon);
	    logp+=logdpoisson(img(i,j),param.bkg+param.ampl*pow(1+r_sq/(param.r0*param.r0),-3*param.beta+.5))+log(expmap(i,j));
	  }
      }
    //cout<<logp<<endl;
    return logp;
  }

  void do_var_range(double& x1,double& x2,const variable& x0,size_t ndim)const
  {
    switch(ndim)
      {
      case 0:
	x1=.001;
	x2=100;
	break;
      case 1:
      case 2:
	x1=0;
	x2=2048;
	break;
      case 3:
	x1=-1;
	x2=1;
	break;
      case 4:
	x1=0;
	x2=2*atan(1)*4;
	break;
      case 5:
	x1=0;
	x2=10;
	break;
      case 6:
	x1=0.001;
	x2=2;
	break;
      case 7:
	x1=0.00;
	x2=1;
	break;
      default:
	assert(0);
	break;
      }
  }

  
};


int main()
{
  prob model("img.fits","exp.fits");
  variable init_var;
  
  
  //cout<<sn.eval(init_var)<<endl;
  
  int s=(time(0));
  
  //cout<<"seed:"<<s<<endl;
  srand(s);
  //return 0;
  variable var;
  u_random<double> rng;
  for(int i=0;i<1000;++i)
    {
      gibbs_sample(model,var,true,rng,1,true);

      cout<<i<<" ";
      //cout<<init_var[1]<<" "<<init_var[2]<<endl;
      for(int j=0;j<var.size();++j)
	{
	  cout<<var[j]<<" ";
	}
      cout<<endl;
    }
}
