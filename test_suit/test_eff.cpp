#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/special_function.hpp>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cassert>

using namespace std;
using namespace mcmc_utilities;

class eff_distribution
  :public probability_density_md<double,std::vector<double> >
{
private:
  std::vector<double> E;
  std::vector<int> nrec;
  std::vector<int> ninj;
public:
  eff_distribution()
  {
    ifstream ifs("eff.txt");
    for(;;)
      {
	double e1;
	int nrec1,ninj1;
	ifs>>e1>>nrec1>>ninj1;
	if(!ifs.good())
	  {
	    break;
	  }
	E.push_back(e1);
	nrec.push_back(nrec1);
	ninj.push_back(ninj1);
      }
  }

  double do_eval_log(const std::vector<double>& x)const
  {
    double A=x.at(0);
    double B=x.at(1);
    double mu=x.at(2);
    double sigma=x.at(3);
    double log_p=0;

    for(unsigned int i=0;i<E.size();++i)
      {
	double eff=A+(B-A)*phi((E[i]-mu)/sigma);
	if(eff>=1)
	  {
	    cout<<"A="<<A<<endl;
	    cout<<"B="<<B<<endl;
	    cout<<"mu="<<mu<<endl;
	    cout<<"sigma="<<sigma<<endl;

	  }
	//dbin<double,double> d(eff,ninj[i]);
	double log_p1=logdbin(nrec[i],eff,ninj[i]);
	log_p+=log_p1;
      }
    return log_p;
  }

  void do_var_range(double& x1,double& x2,const std::vector<double>& x,size_t ndim)const
  {
    if(ndim==0||ndim==1)
      {
	x1=.001;x2=1;
      }
    else if(ndim==2||ndim==3)
      {
	x1=.001;x2=100;
      }
  }

  
  
};

int main()
{
  eff_distribution cd;
  std::vector<double> x;
  x.push_back(.5);
  x.push_back(1);
  x.push_back(15);
  x.push_back(12);
  //cout<<cd.eval_log(x)<<endl;
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,x,1,u_random<double>);
      //if(n>100)
	{
	  for(unsigned int i=0;i<x.size();++i)
	    {
	      cout<<x[i]<<" ";
	    }
	  cout<<endl;
	}
    }
}
