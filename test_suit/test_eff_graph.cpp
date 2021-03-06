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
#include <typeinfo>
using namespace std;
using namespace mcmc_utilities;

template <typename T>
using std_vector=std::vector<T>;

class data_loader
{
public:
  std_vector<double> vec_e,vec_nrec,vec_ninj;
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
  graph<double,std::string,std_vector> g;
  data_loader dl("eff.txt");

  auto pA=std::shared_ptr<node<double,std_vector> >(new fixed_uniform_node<double,std_vector>(.001,1-1e-5));
  auto pB=std::shared_ptr<node<double,std_vector> >(new fixed_uniform_node<double,std_vector>(.001,1-1e-5));
  auto pmu=std::shared_ptr<node<double,std_vector> >(new fixed_uniform_node<double,std_vector>(.001,100-1e-5));
  auto psigma=std::shared_ptr<node<double,std_vector> >(new fixed_uniform_node<double,std_vector>(.001,100-1e-5));
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
      std_vector<std::pair<std::shared_ptr<node<double,std_vector> >,size_t> > pp{{pA,0},{pB,0},{pE,0},{pmu,0},{psigma,0}};
      auto p_eff=new str_node<double,std_vector>("A+(B-A)*phi((E-mu)/sigma)",{"A","B","E","mu","sigma"});
      g.add_node(p_eff,tag_eff,pp);

      
      std::string tag_nrec="nrec";
      tag_nrec+=std::to_string(i);
      g.add_node(dl.get_nrec(i),tag_nrec,{{tag_eff,0},{tag_ninj,0}});
    }

  std::cerr<<"*********"<<std::endl;
  
  graph<double,std::string,std_vector> g2;
  g2.copy_from(g);

  auto topology=g2.topology();

  ofstream ofs("eff_topology.txt");

  for(auto& p:topology)
    {
      std::string tag=p.first;
      std_vector<std::pair<std::string,size_t> > parents=p.second;
      ofs<<tag;
      for(auto q:parents)
	{
	  ofs<<" ("<<q.first<<" , "<<q.second<<" )";
	}
      ofs<<endl;
    }
  ofs.close();
  auto A=g2.get_monitor("A",0);
  auto B=g2.get_monitor("B",0);
  auto mu=g2.get_monitor("mu",0);
  auto sigma=g2.get_monitor("sigma",0);
  g2.set_value("A",0,.01);
  g2.set_value("B",0,.999506);
  g2.set_value("mu",0,13);
  g2.set_value("sigma",0,17);
  g2.initialize();

  
  for(int i=0;i<30000;++i)
    {
      g2.sample(rnd1);
      cout<<A()<<" "<<B()<<" "<<mu()<<" "<<sigma()<<endl;
    }
}
