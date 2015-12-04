#include <core/graph.hpp>
#include <nodes/poisson_node.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/const_node.hpp>
#include <iostream>
#include <core/urand.hpp>

using namespace std;
using namespace mcmc_utilities;

class rnd
  :public base_urand<double>
{
private:
  double do_rand()const
  {
    return rand()/(double)RAND_MAX;
  }
}rnd1;



int main()
{
  auto p_poisson=new poisson_node<double>;
  auto p_lambda=std::shared_ptr<node<double> >(new const_node<double>(100));
  
  graph<double,std::string> g;

  g.add_node(p_lambda,"lambda");
  g.add_node(p_poisson,"poisson",{{p_lambda,0}});
  g.set_value("poisson",0,5);
  auto m=g.get_monitor("poisson",0);

  for(int i=0;i<10000;++i)
    {
      g.sample(rnd1);
      cout<<m()<<endl;
    }
}
