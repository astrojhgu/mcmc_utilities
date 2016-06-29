#include <core/graph.hpp>
#include <nodes/poisson_node.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/const_node.hpp>
#include <iostream>
#include <core/urand.hpp>

using namespace std;
using namespace mcmc_utilities;

template <typename T>
using std_vector=std::vector<T>;

class rnd
  :public base_urand<double>
{
private:
  double do_rand()override
  {
    return rand()/(double)RAND_MAX;
  }
}rnd1;



int main()
{
  auto p_poisson=new poisson_node<double,std_vector>;
  auto p_lambda=std::shared_ptr<node<double,std_vector> >(new const_node<double,std_vector>(100));
  
  graph<double,std::string,std_vector> g;

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
