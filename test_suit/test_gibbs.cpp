#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/distributions.hpp>
#include <random/uniform.h>
#include <core/urand.hpp>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cassert>

using namespace std;
using namespace mcmc_utilities;

class test_distribution
  :public probability_density_md<double,std::vector<double> >
{
public:
  double do_eval_log(const std::vector<double>& x)const
  {
    double log_p=-(x[0]*x[0]+x[1]*x[1]);
    return log_p;
  }

  std::pair<double,double> do_var_range(const std::vector<double>& x0,size_t ndim)const
  {
    return make_pair(-20.,20.);
  }


  std::vector<double> do_candidate_points(const std::vector<double>& x0,size_t ndim)const
  {
    if(ndim==2)
      {
	std::vector<double> res;
	for(double x=-5;x<5;x+=.5)
	  {
	    res.push_back(x);
	  }
	return res;
      }
    return std::vector<double>();
    
  }
};

int main()
{
  test_distribution cd;
  std::vector<double> x(2);
  urand<double> rng;
  gibbs_sample(cd,x,rng);
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,x,rng);
      if(n>100)
	{
	  for(int i=0;i<2;++i)
	    {
	      cout<<x[i]<<" ";
	    }
	  cout<<endl;
	}
    }
}
