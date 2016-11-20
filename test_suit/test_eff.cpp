#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/distributions.hpp>
#include <core/urand.hpp>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cassert>
#include <core/ensemble_sample.hpp>
#include <rng/prng.hpp>
template <typename T>
using std_vector=std::vector<T>;

using namespace std;
using namespace mcmc_utilities;
const double pi=4*atan(1);
struct variable
{
  typedef double value_type;
  double data[4];
  
  double& A;
  double& B;
  double& mu;
  double& sigma;

  int size()const
  {
    return 4;
  }
  
  variable()
    :A(data[0]),B(data[1]),mu(data[2]),sigma(data[3])
  {
    for(int i=0;i<4;++i)
      {
	data[i]=0;
      }
  }

  variable(const variable& rhs)
    :A(data[0]),B(data[1]),mu(data[2]),sigma(data[3])
  {
    for(int i=0;i<4;++i)
      {
	data[i]=rhs.data[i];
      }
  }

  variable& operator=(const variable& rhs)
  {
    for(int i=0;i<4;++i)
      {
	data[i]=rhs.data[i];
      }
    return *this;
  }

  double& operator[](size_t i)
  {
    return data[i];
  }

  const double& operator[](size_t i)const
  {
    return data[i];
  }
};


class eff_distribution
  :public probability_density_md<double,variable,std_vector>
{
private:
  std::vector<double> E;
  std::vector<int> nrec;
  std::vector<int> ninj;
public:
  eff_distribution()
    :E(),nrec(),ninj()
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

  double do_eval_log(const variable& x,int n)const
  {
    if(x.A<0||x.A>1||x.B<0||x.B>1)
      {
	return -std::numeric_limits<double>::infinity();
      }
    
    double log_p=0;

    for(unsigned int i=0;i<E.size();++i)
      {
	double eff=x.A+(x.B-x.A)*phi((E[i]-x.mu)/x.sigma);
	//double eff=x.A+(x.B-x.A)*(atan((E[i]-x.mu)/x.sigma)+(pi/2))/pi;
	//dbin<double,double> d(eff,ninj[i]);
	double log_p1=logdbin(nrec[i],eff,ninj[i]);
	log_p+=log_p1;
      }
    return log_p;
  }

  std::pair<double,double> do_var_range(const variable& x,size_t ndim)const
  {
    if(ndim==0||ndim==1)
      {
	return make_pair(.001,1-1e-5);
      }
    else
      {
	return make_pair(.001,100);
      }
  }

  
  
};

int main()
{
  eff_distribution cd;
  variable x;
  
  x.A=.5;
  x.B=.75;
  x.mu=15;
  x.sigma=12;
  urand<double> rng;
  prng<double> prng;
  //cout<<cd.eval_log(x)<<endl;

  std::vector<variable> ensemble;
  for(int i=0;i<100;++i)
    {
      variable y;
      y.A=.2+rng()*.2-.1;
      y.B=.9+rng()*.2-.1;
      y.mu=15+rng()*.2-.1;
      y.sigma=12+rng()*.2-.1;
      ensemble.push_back(y);
    }

  
  for(int n=0;n<30000;++n)
    {
      //gibbs_sample(cd,x,rng);
      ensemble=ensemble_sample([&cd](const variable& x){return cd.eval_log(x);},ensemble,prng,4);
      int j;
      do{j=rng()*ensemble.size();}while(j==100);
      
      if(n>100)
	{
	  variable x=ensemble[j];
	  for(unsigned int i=0;i<x.size();++i)
	  {
	    cout<<x[i]<<" ";
	  }
	  cout<<endl;
	}
    }
}
