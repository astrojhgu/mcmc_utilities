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
	if(isnan(log_p1))
	  {
	    cout<<"aa"<<log_p1<<endl;
	    cout<<"bb:"<<eff<<" "<<ninj[i]<<" "<<nrec[i]<<endl;
	    cout.flush();
	    assert(0);
	  }
	log_p+=log_p1;
      }
    if(dbg_status)
      {
	cout<<"logp="<<log_p<<endl;
      }
    return log_p;
  }

  eff_distribution* do_clone()const
  {
    return new eff_distribution(*this);
  }

  void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  {
    x1[0]=0.001;x2[0]=1;
    x1[1]=0.001;x2[1]=1;
    x1[2]=0.001;x2[2]=100;
    x1[3]=0.001;x2[3]=100;
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
