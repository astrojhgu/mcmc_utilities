#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
#include <math/functions.hpp>
#include <helper/vnode.hpp>
#include <helper/graph_builder.hpp>
#include <nodes/normal_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/arithmetic_node.hpp>
#include <nodes/obs_bin_node.hpp>
#include <nodes/obs_const_node.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/phi_node.hpp>
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
  

  //gb.add_node(new_node<double,double>("cf")=);
  
  //gb.add_node(new_node<double,double>("cf")=const_vnode<double,double>(.5));
  //  gb.add_node(new_node<double,double>("n")=normal_vnode<double,double>({cpnt<double,double>(const_vnode<double,double>(.6))+cpnt<double,double>(const_vnode<double,double>(.6))+cpnt<double,double>(const_vnode<double,double>(.6))*cpnt<double,double>(const_vnode<double,double>(.6)),0},{const_vnode<double,double>("d",1),0}));

  //gb.add_node(new_node<double,double>("n2")=normal_vnode<double,double>({const_vnode<double,double>("m",0),0},{vnode<double,double>("d"),0}));
  std::vector<double> E_vec{2,6,10,14,18,22,26,30,70,38,42,46,50,54,58,62,66,74,78,82,86,90,94,98,34};
  std::vector<double> nrec_vec{23,71,115,159,200,221,291,244,44,221,210,182,136,119,79,81,61,41,32,32,31,22,18,11,277};
  std::vector<double> ninj_vec{96,239,295,327,345,316,349,281,45,235,217,185,140,121,79,81,61,41,32,32,31,22,18,11,298};


  //gb.add_node("A",uniform_vnode<double,double>("B",0,1)+uniform_vnode<double,double>("C",0,1));
  //auto p=uniform_vnode<double,double>("B",0,1)+uniform_vnode<double,double>("C",0,1);
  //cout<<typeid(*p.parents[0].first.get()).name()<<endl;
  //cout<<typeid(*uniform_vnode<double,double>("C",0,1).clone().get()).name()<<endl;
  
  
  gb.add_node("eff",
	      uniform_vnode<double,double>("A",0,1)(0)+
	      (uniform_vnode<double,double>("B",0,1)
	       -
	       vnode<double,double>("A"))
	      *
	      phi<double,double>(
				 (obs_const_vnode<double,double>("E",E_vec)-uniform_vnode<double,double>("mu",0.,100.))/
				 uniform_vnode<double,double>("sigma",0,100))
	      
	      );

  #if 1
  gb.add_node("nrec",
	      obs_bin_vnode<double,double>(
					   nrec_vec,{vnode<double,double>("eff"),0},{obs_const_vnode<double,double>("ninj",ninj_vec),0}
					   )
	      );
	       
  #endif
  
  gb.validate();
  

  
  ofstream ofs("a.gv");
  graph2gv(gb,ofs);

  //return 0;
  graph<double,double,std::string> g;
  gb.build(g);
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
