#include <core/arms.hpp>
#include <core/urand.hpp>
using namespace mcmc_utilities;

class normal_dist
  :public probability_density_1d<double>
{
public:
  double do_eval_log(const double& x)const
  {
    return -(x-1)*(x-1)/(2*1e-9*1e-9);
  }

  std::pair<double,double> do_var_range()const
  {
    return std::pair<double,double>(-1000,1000);
  }
};
  
urand<double> rng;

int main()
{
  normal_dist nd;
  size_t n=0;
  for(int i=0;i<1000;++i)
    {
      std::cout<<arms(nd,500.0,10,rng,n)-1<<std::endl;
    }
}

