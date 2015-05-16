#include <core/distribution.hpp>
#include <core/gibbs_sampler.hpp>
#include <math/special_function.hpp>
#include <random/uniform.h>
//#include <distribution/dbin.hpp>
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
    double log_p=-(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    return log_p;
  }

  test_distribution* do_clone()const
  {
    return new test_distribution(*this);
  }

  void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  {
    x1.resize(3);
    x2.resize(3);
    x1[0]=x1[1]=x1[2]=-20;
    x2[0]=x2[1]=x2[2]=20;
  }  
};

int main()
{
  test_distribution cd;
  std::vector<double> x(3);
  ranlib::Uniform<double> uf;
  uf.seed(time(0));
  //cout<<cd.eval_log(x)<<endl;
  gibbs_sample(cd,x,1,[&uf](){return uf.random();},2);
  for(int n=0;n<10000;++n)
    {
      gibbs_sample(cd,x,0,[&uf](){return uf.random();},2);
      if(n>100)
	{
	  for(int i=0;i<3;++i)
	    {
	      cout<<x[i]<<" ";
	    }
	  cout<<endl;
	}
    }
}
