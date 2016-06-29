#include <core/graph.hpp>
#include <nodes/poisson_node.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/str_node.hpp>
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
  auto pY=std::shared_ptr<node<double,std_vector> >(new uniform_node<double,std_vector>(-1000,1000));
  auto p_poisson=new poisson_node<double,std_vector>;
  auto pl=std::shared_ptr<node<double,std_vector> >(p_poisson);
  auto pth=std::shared_ptr<node<double,std_vector> >(new const_node<double,std_vector>(150));
  
  auto pc=std::shared_ptr<node<double,std_vector> >(new str_node<double,std_vector>("log(1+(y-theta)^2+1)",{"y","theta"}));

  p_poisson->set_observed(0,true);
  p_poisson->set_value(0,0);
  graph<double,std::string,std_vector> g;
  g.add_node(pY,"Y");
  g.add_node(pth,"theta");
  g.add_node(pc,"cauthy",{{pY,0},{pth,0}});
  g.add_node(pl,"l",{{pc,0}});
  g.initialize();
  auto m=g.get_monitor("Y",0);

  for(int i=0;i<100000;++i)
    {
      g.sample(rnd1);
      cout<<m()<<endl;
    }
}
