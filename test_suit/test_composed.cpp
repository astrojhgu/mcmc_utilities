#include <core/composed_node.hpp>
#include <nodes/arithmetic_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/str_node.hpp>
#include <iostream>

using namespace std;
using namespace mcmc_utilities;

int main()
{
  /*
  composed_node<double,std::string> cn(3,1);
  cn.add_node(new add_node<double>(),"a",{{"",0},{"",0}},{});
  cn.add_node(new add_node<double>(),"d",{{"a",0},{"",0}},{0});
  const_node<double> a1(4);
  const_node<double> a2(3);
  const_node<double> a3(2);
  cn.connect_to_parent(&a1,0,0);
  cn.connect_to_parent(&a2,1,0);
  cn.connect_to_parent(&a3,2,0);
  cout<<cn.value(0)<<endl;;
  */

  //auto p=compose_node<double>("sin(a+a)*A+b",{"a","b","A"});

  str_node<double> nn("a+b-a-b-c",{"a","b","c"});

  str_node<double> nn1("a",{"a"});
  str_node<double> nn2("a",{"a"});
  nn.connect_to_parent(new const_node<double>(4.0),0,0);
  nn.connect_to_parent(new const_node<double>(2.0),1,0);
  nn.connect_to_parent(new const_node<double>(100003.0),2,0);
  nn1.connect_to_parent(&nn,0,0);
  nn2.connect_to_parent(&nn1,0,0);

  cout<<nn2.value(0)<<endl;
}
