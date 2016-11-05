#include <core/mixture_node.hpp>
#include <core/graph.hpp>
#include <core/urand.hpp>
#include <core/tag_t.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/normal_node.hpp>
#include <nodes/bern_node.hpp>
#include <vector>
#include <iostream>
#include <nodes/const_node.hpp>

template <typename T>
using std_vector=std::vector<T>;

using namespace std;
using namespace mcmc_utilities;
urand<double> rnd1;

std::vector<double> generate_data()
{
  graph<double,tag_t,std_vector> g;
  
  auto mn=new mixture_node<double,std_vector>({std::shared_ptr<stochastic_node<double,std_vector> >(new normal_node<double,std_vector>()),std::shared_ptr<stochastic_node<double,std_vector> >(new normal_node<double,std_vector>())});
  auto mu1=new const_node<double,std_vector>(5.2);
  auto sigma1=new const_node<double,std_vector>(1);

  auto mu2=new const_node<double,std_vector>(1);
  auto sigma2=new const_node<double,std_vector>(1);
  //const_node<double,std_vector> id(1);

  auto p=new const_node<double,std_vector>(.2);
  auto id=new bern_node<double,std_vector>;
  //id.connect_to_parent(&p,0,0);
  //mn.connect_to_parent(&id,0,0);
  //mn.connect_to_parent(&mu1,1,0);
  //mn.connect_to_parent(&sigma1,2,0);
  //mn.connect_to_parent(&mu2,3,0);
  //mn.connect_to_parent(&sigma2,4,0);

  g.add_node(p,{"p"});
  g.add_node(id,{"id"},{{tag_t("p"),0}});
  g.add_node(mu1,{"mu1"});
  g.add_node(sigma1,{"sigma1"});
  g.add_node(mu2,{"mu2"});
  g.add_node(sigma2,{"sigma2"});
  g.add_node(mn,{"mn"},{{tag_t("id"),0},{tag_t("mu1"),0},{tag_t("sigma1"),0},{tag_t("mu2"),0},{tag_t("sigma2"),0}});
  auto mon=g.get_monitor(tag_t("mn"),0);
  auto mon2=g.get_monitor(tag_t("id"),0);
  std::vector<double> result;
  for(int i=0;i<1000;++i)
    {
      g.sample(rnd1);
      double x=mon();
      result.push_back(mon());
      //std::cout<<x<<std::endl;
    }
  return result;
}



int main()
{
  srand(time(0));
  
  auto data=generate_data();
  //return 0;
  graph<double,tag_t,std_vector> g;
  auto p=new fixed_uniform_node<double,std_vector>(1e-6,1-1e-6);
  //p->set_observed_value(0,.2);
  //auto mu1=new const_node<double,std_vector>(5.2);
  //auto sigma1=new const_node<double,std_vector>(1);

  //auto mu2=new const_node<double,std_vector>(1);
  //auto sigma2=new const_node<double,std_vector>(.5);
  
  auto mu1=new fixed_uniform_node<double,std_vector>(0,7);
  auto sigma1=new fixed_uniform_node<double,std_vector>(.01,10);

  auto mu2=new fixed_uniform_node<double,std_vector>(0,7);
  auto sigma2=new fixed_uniform_node<double,std_vector>(.01,10);
  g.add_node(p,{"p_prior"});
  g.add_node(mu1,{"mu1"});
  g.add_node(sigma1,{"sigma1"});
  g.add_node(mu2,{"mu2"});
  g.add_node(sigma2,{"sigma2"});
  
  for(int i=0;i<data.size();++i)
    {
      auto id=new bern_node<double,std_vector>;
      g.add_node(id,{"id",i},{{tag_t("p_prior"),0}});
      auto mn=new mixture_node<double,std_vector>({std::shared_ptr<stochastic_node<double,std_vector> >(new normal_node<double,std_vector>()),std::shared_ptr<stochastic_node<double,std_vector> >(new normal_node<double,std_vector>())});
      mn->set_observed_value(0,data[i]);
      g.add_node(mn,{"mn",i},{{tag_t("id",i),0},{tag_t("mu1"),0},{tag_t("sigma1"),0},{tag_t("mu2"),0},{tag_t("sigma2"),0}});
      //std::cerr<<"observed:"<<mn->is_observed(0)<<std::endl;
    }
  g.freeze_topology();
  
  auto mon_id=g.get_monitor(tag_t("p_prior"),0);
  auto mon_mu1=g.get_monitor(tag_t("mu1"),0);
  auto mon_sigma1=g.get_monitor(tag_t("sigma1"),0);
  auto mon_mu2=g.get_monitor(tag_t("mu2"),0);
  auto mon_sigma2=g.get_monitor(tag_t("sigma2"),0);
  std::cout<<"p mu1 sigma1 mu2 sigma2"<<std::endl;
  for(int i=0;i<1000;++i)
    {
      g.sample(rnd1);
      std::cout<<mon_id()<<" "<<mon_mu1()<<" "<<mon_sigma1()<<" "<<mon_mu2()<<" "<<mon_sigma2()<<std::endl;
    }
}
