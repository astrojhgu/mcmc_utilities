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


class data_loader
{
public:
  std::vector<double> vec_e,vec_nrec,vec_ninj;
  data_loader(const char* fname)
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

  shared_ptr<node<double> > get_energy(int i)const
  {
    auto p=new const_node<double>(vec_e[i]);
    //p->set_observed(0,true);
    return std::shared_ptr<node<double> >(p);
  }

  shared_ptr<node<double> > get_ninj(int i)const
  {
    auto p=new const_node<double>(vec_ninj[i]);
    //p->set_observed(0,true);
    return std::shared_ptr<node<double> >(p);
  }

  shared_ptr<node<double> > get_nrec(int i)const
  {
    auto p=new bin_node<double>();
    p->set_value(0,vec_nrec[i]);
    p->set_observed(0,true);
    return shared_ptr<node<double> >(p);
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
  graph<double,std::string> g;
  data_loader dl("eff.txt");

  auto pMu=std::shared_ptr<node<double> >(new const_node<double>(0));
  auto pSigma=std::shared_ptr<node<double> >(new const_node<double>(1));
  auto pk=std::shared_ptr<node<double> >(new const_node<double>(1));
  auto pT1=std::shared_ptr<node<double> >(new t_node<double>());
  auto pT2=std::shared_ptr<node<double> >(new t_node<double>());
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
