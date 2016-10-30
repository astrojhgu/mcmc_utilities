#include <vector>
#include <iostream>
#include <core/graph.hpp>
#include <core/tag_t.hpp>
#include <core/urand.hpp>
#include <nodes/cond_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/normal_node.hpp>
#include <nodes/bern_node.hpp>
#include <core/likelihood_adapter.hpp>
#include <nodes/children_selector_node.hpp>
using namespace std;
using namespace mcmc_utilities;

template <typename T>
using std_vector=std::vector<T>;
int main()
{
  urand<double> rnd1;
  graph<double,tag_t,std_vector> g1;
  g1.add_node(new const_node<double,std_vector>(.01),{"p"});
  g1.add_node(new bern_node<double,std_vector>,tag_t("bern"),{{tag_t("p"),0}});
  g1.add_node(new const_node<double,std_vector>(0),{"m1"});
  g1.add_node(new const_node<double,std_vector>(.2),{"s1"});
  g1.add_node(new normal_node<double,std_vector>,{"g1"},{{tag_t("m1"),0},{tag_t("s1"),0}});
  g1.add_node(new const_node<double,std_vector>(2),{"m2"});
  g1.add_node(new const_node<double,std_vector>(.2),{"s2"});
  g1.add_node(new normal_node<double,std_vector>,{"g2"},{{tag_t("m2"),0},{tag_t("s2"),0}});

  g1.add_node(new cond_node<double,std_vector>,{"result"},{{tag_t("bern"),0},{tag_t("g1"),0},{tag_t("g2"),0}});
  g1.initialize();

  auto m=g1.get_monitor(tag_t("result"),0);
  auto mb=g1.get_monitor(tag_t("bern"),0);
  std::vector<double> simulated_data;
  for(int i=0;i<1000;++i)
    {
      g1.sample(rnd1);
      double x=m();
      simulated_data.push_back(x);
      //std::cout<<x<<std::endl;
    }

  graph<double,tag_t,std_vector> g2;
  /*
  g2.add_node(new const_node<double,std_vector>(0),{"m1"});
  g2.add_node(new const_node<double,std_vector>(.2),{"s1"});
   
  g2.add_node(new const_node<double,std_vector>(2),{"m2"});
  g2.add_node(new const_node<double,std_vector>(.2),{"s2"});
  */  
  auto m1=new uniform_node<double,std_vector>(-10,10);
  auto s1=new uniform_node<double,std_vector>(.1,10);
  auto m2=new uniform_node<double,std_vector>(-10,10);
  auto s2=new uniform_node<double,std_vector>(.1,10);
  m1->set_value(0,0);
  s1->set_value(0,.2);
  m2->set_value(0,2);
  s2->set_value(0,.2);
  g2.add_node(m1,{"m1"});
  g2.add_node(s1,{"s1"});
   
  g2.add_node(m2,{"m2"});
  g2.add_node(s2,{"s2"});


  auto p=new uniform_node<double,std_vector>(1e-5,1-1e-5);
  p->set_value(0,.2);
  g2.add_node(p,{"p"});
  
  for(int i=0;i<simulated_data.size();++i)
    {
      auto v=new const_node<double,std_vector>(simulated_data[i]);
      g2.add_node(new bern_node<double,std_vector>,tag_t("bern",i),{{tag_t("p"),0}});
      auto nn1=new likelihood_adapter<double,std_vector>(new normal_node<double,std_vector>);
      auto nn2=new likelihood_adapter<double,std_vector>(new normal_node<double,std_vector>);
      auto csn=new children_selector_node<double,std_vector>(2);
      g2.add_node(v,tag_t("obsx",i));
      g2.add_node(csn,{"csn",i},{{tag_t("bern",i),0},{tag_t("obsx",i),0}});
      g2.add_node(nn1,{"nn1",i},{{tag_t("m1"),0},{tag_t("s1"),0},{tag_t("csn",i),0}});
      g2.add_node(nn2,{"nn2",i},{{tag_t("m2"),0},{tag_t("s2"),0},{tag_t("csn",i),1}});
    }
  auto mp=g2.get_monitor(tag_t("p"),0);
  
  g2.initialize();
  for(int i=0;i<10000;++i)
    {
      g2.sample(rnd1);
      std::cout<<mp()<<std::endl;
    }
}
