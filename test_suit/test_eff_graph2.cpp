#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
#include <math/functions.hpp>
#include <helper/vnode.hpp>
#include <nodes/unary_node.hpp>
#include <helper/graph_builder.hpp>
#include <nodes/normal_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/arithmetic_node.hpp>
#include <nodes/obs_bin_node.hpp>
#include <nodes/obs_const_node.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/phi_node.hpp>
#include <nodes/bvnormal_node.hpp>
#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <atomic>

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
  graph_builder<double,double> gb;
  

  std::vector<double> E_vec{2,6,10,14,18,22,26,30,70,38,42,46,50,54,58,62,66,74,78,82,86,90,94,98,34};
  std::vector<double> nrec_vec{23,71,115,159,200,221,291,244,44,221,210,182,136,119,79,81,61,41,32,32,31,22,18,11,277};
  std::vector<double> ninj_vec{96,239,295,327,345,316,349,281,45,235,217,185,140,121,79,81,61,41,32,32,31,22,18,11,298};

  gb.add_node(vuniform("A",0,1));
  gb.add_node(vuniform("B",0,1));
  gb.add_node(vuniform("mu",0.,100.));
  gb.add_node(vuniform("sigma",0.,100.));
  gb.add_node("eff",
	      //vn("A")+(vn("B")-vn("A"))*vphi((obs_vconst("E",E_vec)-vn("mu"))/vn("sigma"))
	      vn("A")+(vn("B")-vn("A"))*vunary((obs_vconst("E",E_vec)-vn("mu"))/vn("sigma"),{[](const double& x){return phi(x);}})
	      );
  gb.add_node("ninj",obs_vconst("ninj",ninj_vec));
  gb.add_node("nrec",
	      obs_vbin(
		       nrec_vec,{vn("eff"),0},{vn("ninj"),0}
		       )
	      );
  gb.validate();
  
  ofstream ofs("a.gv");
  graph2dot1(gb,ofs);


  graph<double,double,std::string> g;
  gb.build(g);
  g.set_value("A",0,5.14989e-05);
  g.set_value("B",0,0.999406);
  g.set_value("mu",0,13.2949);
  //g.set_value("mu",0,31.7934);
  //g.set_value("sigma",0,17.5035);
  g.set_value("sigma",0,17.6002);

  
  g.sample(rnd1);
  //return 0;
  for(int i=0;i<10;++i)
    {
      g.sample(rnd1);
    }
  
  for(int i=0;i<30000;++i)
    {
      g.sample(rnd1);
      auto p=g.get_params();
      for(auto j : p)
	{
	  cout<<j<<" ";
	}
      cout<<endl;
      //  cout<<g.get_value("n",0)<<" "<<g.get_value("n2",0)<<endl;
    }
}
