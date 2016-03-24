#include <core/arms.hpp>
#include <core/urand.hpp>
using namespace mcmc_utilities;


urand<double> rng;

int main()
{
  
  size_t n=0;
  for(int i=0;i<1000;++i)
    {
      std::cout<<arms([](const double& x){return -(x-1)*(x-1)/(2*1e-9*1e-9);},{-1000.,1000.},{-100.,-.001,0,.001,100.},500.0,10,rng,n)-1<<std::endl;
    }
}

