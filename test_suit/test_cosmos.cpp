#include <cmath>
#include <fstream>
#include <iostream>
#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
#include <core/urand.hpp>
#include <core/tag_t.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/bin_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/str_node.hpp>
#include <nodes/normal_node.hpp>
#include <tools/dump_graph_topology.hpp>
#include <math/functions.hpp>
#include <io/jags_data_reader.hpp>
#include <core/urand.hpp>

using namespace std;
using namespace mcmc_utilities;

template <typename T>
using std_vector=std::vector<T>;


int main()
{
  urand<double> rnd;
  graph<double,tag_t,std_vector> g;
  ifstream ifs("data8.12.1.dat.R");
  auto data=read_jags_data<double>(ifs);

  auto obsz=data.find("obsz")->second;
  auto errz=data.find("errz")->second;
  auto obsm=data.find("obsm")->second;
  auto errmag=data.find("errmag")->second;
  auto obsc=data.find("obsc")->second;
  auto errobsc=data.find("errobsc")->second;
  auto obsx=data.find("obsx")->second;
  auto errobsx=data.find("errobsx")->second;

  std::shared_ptr<node<double,std_vector> > pMm(new fixed_uniform_node<double,std_vector>(-20.3,-18.3));
  std::shared_ptr<node<double,std_vector> > palpha(new fixed_uniform_node<double,std_vector>(-2,2));
  std::shared_ptr<node<double,std_vector> > pbeta(new fixed_uniform_node<double,std_vector>(-4,4));
  std::shared_ptr<node<double,std_vector> > pcm(new fixed_uniform_node<double,std_vector>(-3,3));
  std::shared_ptr<node<double,std_vector> > pxm(new fixed_uniform_node<double,std_vector>(-10,10));
  std::shared_ptr<node<double,std_vector> > plgintrscatM(new fixed_uniform_node<double,std_vector>(-3,0));
  std::shared_ptr<node<double,std_vector> > pintrscatM(new str_node<double,std_vector>("10^x",{"x"}));
  std::shared_ptr<node<double,std_vector> > plgintrscatx(new fixed_uniform_node<double,std_vector>(-5,2));
  std::shared_ptr<node<double,std_vector> > pintrscatx(new str_node<double,std_vector>("10^x",{"x"}));
  std::shared_ptr<node<double,std_vector> > plgintrscatC(new fixed_uniform_node<double,std_vector>(-5,2));
  std::shared_ptr<node<double,std_vector> > pintrscatC(new str_node<double,std_vector>("10^x",{"x"}));
  std::shared_ptr<node<double,std_vector> > pH0_mu(new const_node<double,std_vector>(72));
  std::shared_ptr<node<double,std_vector> > pH0_sigma(new const_node<double,std_vector>(8.0));
  std::shared_ptr<node<double,std_vector> > pH0(new normal_node<double,std_vector>());
  std::shared_ptr<node<double,std_vector> > pOmega_l(new str_node<double,std_vector>("1-x",{"x"}));
  std::shared_ptr<node<double,std_vector> > pOmega_m(new fixed_uniform_node<double,std_vector>(0,1));
  std::shared_ptr<node<double,std_vector> > pw(new fixed_uniform_node<double,std_vector>(-4,0));

  g.add_node(pMm,{"Mm"});
  g.add_node(palpha,{"alpha"});
  g.add_node(pbeta,{"beta"});
  g.add_node(pcm,{"cm"});
  g.add_node(pxm,{"xm"});
  g.add_node(plgintrscatM,{"lgintrscatM"});
  g.add_node(pintrscatM,{"intrscatM"},{{plgintrscatM,0}});
  
  g.add_node(plgintrscatx,{"lgintrscatx"});
  g.add_node(pintrscatx,{"intrscatx"},{{plgintrscatx,0}});
  
  g.add_node(plgintrscatC,{"lgintrscatC"});
  g.add_node(pintrscatC,{"intrscatC"},{{plgintrscatC,0}});
  g.add_node(pH0_mu,{"H0_mu"});
  g.add_node(pH0_sigma,{"H0_sigma"});
  g.add_node(pH0,{"H0"},{{pH0_mu,0},{pH0_sigma,0}});

  g.add_node(pOmega_m,{"Omega_m"});
  g.add_node(pOmega_l,{"Omega_l"},{{pOmega_m,0}});

  g.add_node(pw,{"w"});

  for(int i=0;i<obsz.size();++i)
    {
      std::shared_ptr<node<double,std_vector> > pz(new fixed_uniform_node<double,std_vector>(1e-6,2));
      g.add_node(pz,{"z",i});
      g.set_value({"z",i},0,obsz[i]);
      std::shared_ptr<node<double,std_vector> > perrz(new const_node<double,std_vector>(errz[i]));
      g.add_node(perrz,{"errz",i});
      auto pobsz=new normal_node<double,std_vector>();
      pobsz->set_value(0,obsz[i]);
      pobsz->set_observed(0,true);
      g.add_node(pobsz,{"obsz",i},{{pz,0},{perrz,0}});

      auto pdistmod=std::shared_ptr<node<double,std_vector> >(new str_node<double,std_vector>("5*log10(D_L(z,H0,Omega_m,1-Omega_m,0,0,w)/3.0857E16)-5",{"z","H0","Omega_m","w"}));

      g.add_node(pdistmod,{"distmod",i},{{pz,0},{pH0,0},{pOmega_m,0},{pw,0}});

      auto px=shared_ptr<node<double,std_vector> >(new normal_node<double,std_vector>());
      g.add_node(px,{"x",i},{{pxm,0},{pintrscatx,0}});
      auto perrobsx=std::shared_ptr<node<double,std_vector> >(new const_node<double,std_vector>(errobsx[i]));
      g.add_node(perrobsx,{"errobsx",i});
      auto pobsx=new normal_node<double,std_vector>();
      pobsx->set_value(0,obsx[i]);
      pobsx->set_observed(0,true);
      g.add_node(pobsx,{"obsx",i},{{px,0},{perrobsx,0}});

      shared_ptr<node<double,std_vector> > pc(new normal_node<double,std_vector>);
      g.add_node(pc,{"c",i},{{pcm,0},{pintrscatC,0}});

      shared_ptr<node<double,std_vector> > perrobs(new const_node<double,std_vector>(errobsc[i]));
      g.add_node(perrobs,{"errobs",i});
      auto pobsc=new normal_node<double,std_vector>;
      pobsc->set_value(0,obsc[i]);
      pobsc->set_observed(0,true);
      g.add_node(pobsc,{"obsc",i},{{pc,0},{perrobs,0}});

      shared_ptr<node<double,std_vector> > pm_mean(new str_node<double,std_vector>("Mm+distmod-alpha*x+beta*c",{"Mm","distmod","alpha","x","beta","c"}));
      g.add_node(pm_mean,{"pm_mean",i},{{pMm,0},{pdistmod,0},{palpha,0},{px,0},{pbeta,0},{pc,0}});

      std::shared_ptr<node<double,std_vector> > pm(new normal_node<double,std_vector>);
      g.add_node(pm,{"m",i},{{pm_mean,0},{pintrscatM,0}});

      std::shared_ptr<node<double,std_vector> > perrmag(new const_node<double,std_vector>(errmag[i]));

      g.add_node(perrmag,{"errmag",i});

      auto pobsm=new normal_node<double,std_vector>;
      pobsm->set_value(0,obsm[i]);
      pobsm->set_observed(0,true);
      g.add_node(pobsm,{"obsm",i},{{pm,0},{perrmag,0}});
    }

  g.initialize();
  auto m_omega_m=g.get_monitor({"Omega_m"},0);
  auto m_w=g.get_monitor({"w"},0);
  auto m_mag=g.get_monitor({"m",23},0);
  auto m_distmod=g.get_monitor({"distmod",23},0);
  auto m_Mm=g.get_monitor({"Mm"},0);
  auto m_z=g.get_monitor({"z",23},0);
  auto m_H0=g.get_monitor({"H0"},0);
  g.set_value({"Omega_m"},0,0.27);
  g.set_value({"w"},0,-1);
  g.set_value({"H0"},0,72);

  ofstream ofs_tp("cosmology.dot");
  topology_dumper<double,std_vector>(g).to_dot(ofs_tp);
  ofs_tp.close();

  auto pdistmod=std::shared_ptr<node<double,std_vector> >(new str_node<double,std_vector>("5*log10(D_L(z,H0,Omega_m,1-Omega_m,0,0,w)/3.0858E16)-5",{"z","H0","Omega_m","w"}));
  //auto pdistmod=std::shared_ptr<node<double,std_vector> >(new str_node<double,std_vector>("D_L(z,H0,Omega_m,1-Omega_m,0,0,w)",{"z","H0","Omega_m","w"}));
  pdistmod->connect_to_parent(new const_node<double,std_vector>(.5),0,0);
  pdistmod->connect_to_parent(new const_node<double,std_vector>(71),1,0);
  pdistmod->connect_to_parent(new const_node<double,std_vector>(0.27),2,0);
  pdistmod->connect_to_parent(new const_node<double,std_vector>(-1),3,0);

  //cout<<pdistmod->value(0)<<endl;
  //return 0;
  //cout<<m_distmod()<<endl;
  //return 0;
  for(int i=0;i<10000;++i)
    {
      g.sample(rnd);
      if(i<100)
	{
	  cerr<<m_omega_m()<<" "<<m_w()<<" "<<m_Mm()<<endl;
	}
      else
	{
	  //cout<<m_omega_m()<<" "<<m_w()<<" "<<m_Mm()<<" "<<m_z()<<" "<<m_H0()<<" "<<m_distmod()<<endl;
	  cout<<m_omega_m()<<" "<<m_w()<<endl;
	}
      //cout<<m_mag()<<endl;
      //cout<<(m_mag()-m_Mm())<<" "<<m_distmod()<<" "<<m_z()<<" "<<m_omega_m()<<" "<<m_w()<<" "<<m_H0()<<endl;
    }
}
