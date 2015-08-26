#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
#include <core/urand.hpp>
#include <math/functions.hpp>
#include <nodes/obs_normal_node.hpp>
#include <nodes/uniform_node.hpp>
#include <cassert>
#include <fstream>
#include <iostream>


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


int main(int argc,char* argv[])
{

  if(argc!=2)
    {
      cerr<<"Usage:"<<argv[0]<<" <data file>"<<endl;
      return -1;
    }
  graph<double,double,std::string> g;
  uniform_node_factory<double,double> unf;
  obs_normal_node_factory<double,double> nmf;

  auto pm1=unf.get_node({0.0,5.0},{});
  auto pm2=unf.get_node({1.0,6.0},{});

  auto pu1=dynamic_pointer_cast<uniform_node<double,double> >(pm1);
  auto pu2=dynamic_pointer_cast<uniform_node<double,double> >(pm2);
  g.add_node(pm1,"A",{});
  g.add_node(pm2,"B",{});

  std::vector<double> data;

  ifstream ifs(argv[1]);

  for(;;)
    {
      double x;
      ifs>>x;
      if(!ifs.good())
	{
	  break;
	}
      data.push_back(x);
    }
  auto pmodel=nmf.get_node({},{data});

  auto pnorm=dynamic_pointer_cast<obs_normal_node<double,double> >(pmodel);
  
  g.add_node(pmodel,"gauss",{{"A",0},{"B",0}});
  pu2->set_value(0,3);
  pu1->set_value(0,2);
  for(double x=1;x<5;x+=.01)
    {
      pu2->set_value(0,x);
      cout<<x<<" "<<pnorm->log_prior_prob()<<endl;
    }



  return 0;
  for(int i=0;i<100;++i)
    {
      g.sample(rnd1);
    }
  

  auto A=g.get_monitor("A",0);
  auto B=g.get_monitor("B",0);
  for(int i=0;i<30000;++i)
    {
      g.sample(rnd1);

      cout<<A()<<" "<<B()<<endl;
      
    }
}
