#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
#include <core/urand.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/bin_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/str_node.hpp>
#include <nodes/t_node.hpp>
#include <math/functions.hpp>
#include <cassert>
#include <fstream>
#include <iostream>

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
  graph<double,std::string,std_vector> g;
  auto pMu=std::shared_ptr<node<double,std_vector> >(new const_node<double,std_vector>(0));
  auto pSigma=std::shared_ptr<node<double,std_vector> >(new const_node<double,std_vector>(1));
  auto pk=std::shared_ptr<node<double,std_vector> >(new const_node<double,std_vector>(1));
  auto pT1=std::shared_ptr<node<double,std_vector> >(new t_node<double,std_vector>());
  auto pT2=std::shared_ptr<node<double,std_vector> >(new t_node<double,std_vector>());
  g.add_node(pMu,"mu");
  g.add_node(pSigma,"sigma");
  g.add_node(pk,"k");
  g.add_node(pT1,"T1",{{pMu,0},{pSigma,0},{pk,0}});
  g.add_node(pT2,"T2",{{pMu,0},{pSigma,0},{pk,0}});
  

  auto T1=g.get_monitor("T1",0);
  auto T2=g.get_monitor("T2",0);
  g.set_value("T1",0,0.);
  g.set_value("T2",0,0.);
  for(int i=0;i<300000;++i)
    {
      g.sample(rnd1);
      rnd1();
      cout<<T1()<<" "<<T2()<<endl;
    }
}
