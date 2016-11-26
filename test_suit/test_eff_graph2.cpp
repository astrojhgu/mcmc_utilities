#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
#include <core/urand.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/bin_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/func_node.hpp>
#include <math/functions.hpp>
#include <core/tag_t.hpp>
#include <tools/dump_graph_topology.hpp>
#include <core/ensemble_sample.hpp>
#include <core/graph_ensemble_sample.hpp>
#include <rng/prng.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <typeinfo>
using namespace std;
using namespace mcmc_utilities;

template <typename T>
using std_vector=std::vector<T>;

class data_loader
{
public:
  std::vector<double> vec_e,vec_nrec,vec_ninj;
  data_loader(const char* fname)
    :vec_e(),
     vec_nrec(),
     vec_ninj()
  {
    ifstream ifs(fname);
    for(;;)
      {
	double e,nrec,ninj;
	ifs>>e>>nrec>>ninj;
	if(!ifs.good())
	  {
	    break;
	  }
	vec_e.push_back(e);
	vec_nrec.push_back(nrec);
	vec_ninj.push_back(ninj);
      }
  }

  shared_ptr<node<double,std_vector> > get_energy(int i)const
  {
    auto p=new const_node<double,std_vector>(vec_e[i]);
    //p->set_observed(0,true);
    return std::shared_ptr<node<double,std_vector> >(p);
  }

  shared_ptr<node<double,std_vector> > get_ninj(int i)const
  {
    auto p=new const_node<double,std_vector>(vec_ninj[i]);
    //p->set_observed(0,true);
    return std::shared_ptr<node<double,std_vector> >(p);
  }

  shared_ptr<node<double,std_vector> > get_nrec(int i)const
  {
    auto p=new bin_node<double,std_vector>();
    p->set_value(0,vec_nrec[i]);
    p->set_observed(0,true);
    return shared_ptr<node<double,std_vector> >(p);
  }

  size_t size()const
  {
    return vec_e.size();
  }
};

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
  prng<double> prng;
  
  graph<double,tag_t,std_vector> g;
  data_loader dl("eff.txt");

  auto pA=std::shared_ptr<node<double,std_vector> >(new fixed_uniform_node<double,std_vector>(.001,1-1e-5));
  auto pB=std::shared_ptr<node<double,std_vector> >(new fixed_uniform_node<double,std_vector>(.001,1-1e-5));
  auto pmu=std::shared_ptr<node<double,std_vector> >(new fixed_uniform_node<double,std_vector>(.001,100-1e-5));
  auto psigma=std::shared_ptr<node<double,std_vector> >(new fixed_uniform_node<double,std_vector>(.001,100-1e-5));

  g.add_node(pA,{"A"});
  g.add_node(pB,{"B"});
  g.add_node(pmu,{"mu"});
  g.add_node(psigma,{"sigma"});
  
  for(int i=0;i<dl.size();++i)
    {
      auto pE=dl.get_energy(i);
      g.add_node(pE,{"E",i});
      
      g.add_node(dl.get_ninj(i),{"ninj",i});

      
      std::vector<std::pair<std::shared_ptr<node<double,std_vector> >,size_t> > pp{{pA,0},{pB,0},{pE,0},{pmu,0},{psigma,0}};
      auto p_eff=new func_node<double,std_vector>([](const std::vector<double>& p)->double
	{
	  double A=p[0];
	  double B=p[1];
	  double E=p[2];
	  double mu=p[3];
	  double sigma=p[4];
	  return A+(B-A)*phi((E-mu)/sigma);
	},5);
      g.add_node(p_eff,{"eff",i},pp);

      g.add_node(dl.get_nrec(i),{"nrec",i},{{tag_t("eff",i),0},{tag_t("ninj",i),0}});
    }

  std::cerr<<"*********"<<std::endl;
  
  graph<double,tag_t,std_vector> g2;
  g2.copy_from(g);

  auto A=g2.get_monitor({"A"},0);
  auto B=g2.get_monitor({"B"},0);
  auto mu=g2.get_monitor({"mu"},0);
  auto sigma=g2.get_monitor({"sigma"},0);
  g2.set_value({"A"},0,.01);
  g2.set_value({"B"},0,.999506);
  g2.set_value({"mu"},0,13);
  g2.set_value({"sigma"},0,17);
  g2.initialize();

  ofstream ofs("eff_topology.dot");
  topology_dumper<double,std_vector>(g2).to_dot(ofs);
  ofs.close();

  std_vector<std_vector<double> > ensemble;
  
  
  for(int i=0;i<30000;++i)
    {
      //g2.sample(rnd1);
      //cout<<A()<<" "<<B()<<" "<<mu()<<" "<<sigma()<<endl;
      ensemble=ensemble_sample(g,ensemble,prng);
      auto p=ensemble[int(urng<double>(rnd1)*8)%8];
      //auto p=ensemble[0];
      if(i<1000)
	{
	  continue;
	}
      for(auto& j:p)
	{
	  cout<<j<<" ";
	}
      cout<<endl;
    }
}
