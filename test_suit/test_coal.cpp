#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <vector>
#include <fstream>
#include <cassert>
using namespace std;
using namespace mcmc_utilities;

int factorial(int n)
{
  assert(n<100);
  if(n==0)
    {
      return 1;
    }
  else
    {
      return n*factorial(n-1);
    }
}
class coal_distribution
  :public probability_density_md<double,std::vector<double> >
{
private:
  std::vector<double> year;
  std::vector<double> n_accident;
public:
  coal_distribution()
  {
    ifstream ifs("coal.data");
    for(;;)
      {
	double y,n;
	ifs>>y>>n;
	if(!ifs.good())
	  {
	    break;
	  }
	year.push_back(y);
	n_accident.push_back(n);
      }
  }

  double do_eval_log(const std::vector<double>& x)const
  {
    double lambda1=x[0];
    double lambda2=x[1];
    double k=x[2];
    double logp=0;
    double lambda=0;
    for(int i=0;i<year.size();++i)
      {
	if(year[i]<k)
	  {
	    lambda=lambda1;
	  }
	else
	  {
	    lambda=lambda2;
	  }
	int n=n_accident[i];
	double logp1=n*log(lambda)-log(factorial(n))-lambda;
	logp+=logp1;
      }
    return logp;
  }

  coal_distribution* do_clone()const
  {
    return new coal_distribution(*this);
  }

  void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  {
    x1[0]=0;x2[0]=100;
    x1[1]=0;x2[1]=100;
    x1[2]=0;x2[2]=year.back()+1;
  }

  
  
};

int main()
{
  coal_distribution cd;
  std::vector<double> x;
  x.push_back(3);
  x.push_back(4);
  x.push_back(50);
  //cout<<cd.eval_log(x)<<endl;
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,x,1,u_random<double>);
      for(int i=0;i<x.size();++i)
	{
	  cout<<x[i]<<" ";
	}
      cout<<endl;
    }
}
