#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
#include <core/urand.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/bin_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/str_node.hpp>
#include <math/functions.hpp>
#include <cassert>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mcmc_utilities;

class eff
  :public deterministic_node<double,double>
{
public:
  eff()
    :deterministic_node<double,double>(5)
  {}

  double do_value(size_t idx)const override
  {
    double A=this->parent(0);
    double B=this->parent(1);
    double E=this->parent(2);
    double mu=this->parent(3);
    double sigma=this->parent(4);
    return A+(B-A)*phi((E-mu)/sigma);
    //return this->parents[0]
  }
};


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

  shared_ptr<node<double,double> > get_energy(int i)const
  {
    auto p=new const_node<double,double>(vec_e[i]);
    //p->set_observed(0,true);
    return std::shared_ptr<node<double,double> >(p);
  }

  shared_ptr<node<double,double> > get_ninj(int i)const
  {
    auto p=new const_node<double,double>(vec_ninj[i]);
    //p->set_observed(0,true);
    return std::shared_ptr<node<double,double> >(p);
  }

  shared_ptr<node<double,double> > get_nrec(int i)const
  {
    auto p=new bin_node<double,double>();
    p->set_value(0,vec_nrec[i]);
    p->set_observed(0,true);
    return shared_ptr<node<double,double> >(p);
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
  double do_rand()const
  {
    return rand()/(double)RAND_MAX;
  }
}rnd1;


int main()
{
  graph<double,double,std::string> g;
  data_loader dl("eff.txt");

  auto pA=std::shared_ptr<node<double,double> >(new uniform_node<double,double>(.001,1-1e-5));
  auto pB=std::shared_ptr<node<double,double> >(new uniform_node<double,double>(.001,1-1e-5));
  auto pmu=std::shared_ptr<node<double,double> >(new uniform_node<double,double>(.001,100-1e-5));
  auto psigma=std::shared_ptr<node<double,double> >(new uniform_node<double,double>(.001,100-1e-5));
  g.add_node(pA,"A");
  g.add_node(pB,"B");
  g.add_node(pmu,"mu");
  g.add_node(psigma,"sigma");
  
  for(int i=0;i<dl.size();++i)
    {
      std::string tag_E("E");
      tag_E+=std::to_string(i);

      auto pE=dl.get_energy(i);
      g.add_node(pE,tag_E);
      
      std::string tag_ninj="ninj";
      tag_ninj+=std::to_string(i);
      g.add_node(dl.get_ninj(i),tag_ninj);

      std::string tag_eff="eff";
      tag_eff+=std::to_string(i);

      //g.add_node(new eff(),tag_eff,{{"A",0},{"B",0},{tag_E,0},{"mu",0},{"sigma",0}});
      std::vector<std::pair<std::shared_ptr<node<double,double> >,size_t> > pp{{pA,0},{pB,0},{pE,0},{pmu,0},{psigma,0}};
      g.add_node(new str_node<double,double>("A+(B-A)*phi((E-mu)/sigma*1)",{"A","B","E","mu","sigma"}),tag_eff,pp);
      std::string tag_nrec="nrec";
      tag_nrec+=std::to_string(i);
      g.add_node(dl.get_nrec(i),tag_nrec,{{tag_eff,0},{tag_ninj,0}});
    }
  

  auto A=g.get_monitor("A",0);
  auto B=g.get_monitor("B",0);
  auto mu=g.get_monitor("mu",0);
  auto sigma=g.get_monitor("sigma",0);
  g.set_value("A",0,.01);
  g.set_value("B",0,.999506);
  g.set_value("mu",0,13);
  g.set_value("sigma",0,17);
  for(int i=0;i<30000;++i)
    {
      g.sample(rnd1);
      cout<<A()<<" "<<B()<<" "<<mu()<<" "<<sigma()<<endl;
    }
}
