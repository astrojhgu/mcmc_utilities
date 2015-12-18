#include <cmath>
#include <fstream>
#include <iostream>
#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
#include <core/urand.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/bin_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/str_node.hpp>
#include <nodes/normal_node.hpp>
#include <math/functions.hpp>
#include <io/jags_data_reader.hpp>
#include <core/urand.hpp>

using namespace std;
using namespace mcmc_utilities;


int main()
{
  urand<double> rnd;
  graph<double,std::string> g;
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

  std::shared_ptr<node<double> > pMm(new uniform_node<double>(-20.3,-18.3));
  std::shared_ptr<node<double> > palpha(new uniform_node<double>(-2,2));
  std::shared_ptr<node<double> > pbeta(new uniform_node<double>(-4,4));
  std::shared_ptr<node<double> > pcm(new uniform_node<double>(-3,3));
  std::shared_ptr<node<double> > pxm(new uniform_node<double>(-10,10));
  std::shared_ptr<node<double> > plgintrscatM(new uniform_node<double>(-3,0));
  std::shared_ptr<node<double> > pintrscatM(new str_node<double>("10^x",{"x"}));
  std::shared_ptr<node<double> > plgintrscatx(new uniform_node<double>(-5,2));
  std::shared_ptr<node<double> > pintrscatx(new str_node<double>("10^x",{"x"}));
  std::shared_ptr<node<double> > plgintrscatC(new uniform_node<double>(-5,2));
  std::shared_ptr<node<double> > pintrscatC(new str_node<double>("10^x",{"x"}));
  std::shared_ptr<node<double> > pH0_mu(new const_node<double>(72));
  std::shared_ptr<node<double> > pH0_sigma(new const_node<double>(8.0));
  std::shared_ptr<node<double> > pH0(new normal_node<double>());
  std::shared_ptr<node<double> > pOmega_l(new str_node<double>("1-x",{"x"}));
  std::shared_ptr<node<double> > pOmega_m(new uniform_node<double>(0,1));
  std::shared_ptr<node<double> > pw(new uniform_node<double>(-4,0));

  g.add_node(pMm,"Mm");
  g.add_node(palpha,"alpha");
  g.add_node(pbeta,"beta");
  g.add_node(pcm,"cm");
  g.add_node(pxm,"xm");
  g.add_node(plgintrscatM,"lgintrscatM");
  g.add_node(pintrscatM,"intrscatM",{{plgintrscatM,0}});
  
  g.add_node(plgintrscatx,"lgintrscatx");
  g.add_node(pintrscatx,"intrscatx",{{plgintrscatx,0}});
  
  g.add_node(plgintrscatC,"lgintrscatC");
  g.add_node(pintrscatC,"intrscatC",{{plgintrscatC,0}});
  g.add_node(pH0_mu,"H0_mu");
  g.add_node(pH0_sigma,"H0_sigma");
  g.add_node(pH0,"H0",{{pH0_mu,0},{pH0_sigma,0}});

  g.add_node(pOmega_m,"Omega_m");
  g.add_node(pOmega_l,"Omega_l",{{pOmega_m,0}});

  g.add_node(pw,"w");

  for(int i=0;i<obsz.size();++i)
    {
      std::shared_ptr<node<double> > pz(new uniform_node<double>(0,2));
      g.add_node(pz,std::string("z")+to_string(i));
      g.set_value(std::string("z")+to_string(i),0,obsz[i]);
      std::shared_ptr<node<double> > perrz(new const_node<double>(errz[i]));
      g.add_node(perrz,std::string("errz")+to_string(i));
      auto pobsz=new normal_node<double>();
      pobsz->set_value(0,obsz[i]);
      pobsz->set_observed(0,true);
      g.add_node(pobsz,std::string("obsz")+to_string(i),{{pz,0},{perrz,0}});

      auto pdistmod=std::shared_ptr<node<double> >(new str_node<double>("5*log10(D_L(z,H0,Omega_m,1-Omega_m,0,0,w)/3.0857E16)-5",{"z","H0","Omega_m","w"}));

      g.add_node(pdistmod,std::string("distmod")+to_string(i),{{pz,0},{pH0,0},{pOmega_m,0},{pw,0}});

      auto px=shared_ptr<node<double> >(new normal_node<double>());
      g.add_node(px,std::string("x")+to_string(i),{{pxm,0},{pintrscatx,0}});
      auto perrobsx=std::shared_ptr<node<double> >(new const_node<double>(errobsx[i]));
      g.add_node(perrobsx,std::string("errobsx")+to_string(i));
      auto pobsx=new normal_node<double>();
      pobsx->set_value(0,obsx[i]);
      pobsx->set_observed(0,true);
      g.add_node(pobsx,string("obsx")+to_string(i),{{px,0},{perrobsx,0}});

      shared_ptr<node<double> > pc(new normal_node<double>);
      g.add_node(pc,std::string("c")+to_string(i),{{pcm,0},{pintrscatC,0}});

      shared_ptr<node<double> > perrobs(new const_node<double>(errobsc[i]));
      g.add_node(perrobs,std::string("errobs")+to_string(i));
      auto pobsc=new normal_node<double>;
      pobsc->set_value(0,obsc[i]);
      pobsc->set_observed(0,true);
      g.add_node(pobsc,std::string("obsc")+to_string(i),{{pc,0},{perrobs,0}});

      shared_ptr<node<double> > pm_mean(new str_node<double>("Mm+distmod-alpha*x+beta*c",{"Mm","distmod","alpha","x","beta","c"}));
      g.add_node(pm_mean,std::string("pm_mean")+to_string(i),{{pMm,0},{pdistmod,0},{palpha,0},{px,0},{pbeta,0},{pc,0}});

      std::shared_ptr<node<double> > pm(new normal_node<double>);
      g.add_node(pm,"m"+to_string(i),{{pm_mean,0},{pintrscatM,0}});

      std::shared_ptr<node<double> > perrmag(new const_node<double>(errmag[i]));

      g.add_node(perrmag,"errmag"+to_string(i));

      auto pobsm=new normal_node<double>;
      pobsm->set_value(0,obsm[i]);
      pobsm->set_observed(0,true);
      g.add_node(pobsm,"obsm"+to_string(i),{{pm,0},{perrmag,0}});
    }

  g.initialize();
  auto m_omega_m=g.get_monitor("Omega_m",0);
  auto m_w=g.get_monitor("w",0);
  auto m_mag=g.get_monitor("m23",0);
  auto m_distmod=g.get_monitor("distmod23",0);
  auto m_Mm=g.get_monitor("Mm",0);
  auto m_z=g.get_monitor("z23",0);
  auto m_H0=g.get_monitor("H0",0);
  g.set_value("Omega_m",0,0.27);
  g.set_value("w",0,-1);
  g.set_value("H0",0,72);


  auto pdistmod=std::shared_ptr<node<double> >(new str_node<double>("5*log10(D_L(z,H0,Omega_m,1-Omega_m,0,0,w)/3.0858E16)-5",{"z","H0","Omega_m","w"}));
  //auto pdistmod=std::shared_ptr<node<double> >(new str_node<double>("D_L(z,H0,Omega_m,1-Omega_m,0,0,w)",{"z","H0","Omega_m","w"}));
  pdistmod->connect_to_parent(new const_node<double>(.5),0,0);
  pdistmod->connect_to_parent(new const_node<double>(71),1,0);
  pdistmod->connect_to_parent(new const_node<double>(0.27),2,0);
  pdistmod->connect_to_parent(new const_node<double>(-1),3,0);

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
	  cout<<m_omega_m()<<" "<<m_w()<<" "<<m_Mm()<<" "<<m_z()<<" "<<m_H0()<<" "<<m_distmod()<<endl;
	}
      //cout<<m_mag()<<endl;
      //cout<<(m_mag()-m_Mm())<<" "<<m_distmod()<<" "<<m_z()<<" "<<m_omega_m()<<" "<<m_w()<<" "<<m_H0()<<endl;
    }
}
