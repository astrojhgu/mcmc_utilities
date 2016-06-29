#include <core/likelihood_adapter.hpp>
#include <nodes/t_node.hpp>
#include <core/graph.hpp>
#include <nodes/uniform_node.hpp>
#include <nodes/const_node.hpp>
#include <core/urand.hpp>
#include <iomanip>
#include <memory>

using namespace mcmc_utilities;
using namespace std;

template <typename T>
using std_vector=std::vector<T>;

int main()
{
  urand<double> rnd;
  graph<double,std::string,std_vector> g;
  auto x_prior=std::shared_ptr<node<double,std_vector> >(new uniform_node<double,std_vector>(-1000,1000));
  auto x_likelihood=std::shared_ptr<node<double,std_vector> >(new likelihood_adapter<double,std_vector>(new t_node<double,std_vector>));
  auto t_mean=std::shared_ptr<node<double,std_vector> >(new const_node<double,std_vector>(1.0));
  auto t_sigma=std::shared_ptr<node<double,std_vector> >(new const_node<double,std_vector>(1.0e-7));
  auto t_dof=std::shared_ptr<node<double,std_vector> >(new const_node<double,std_vector>(1.0));
  g.add_node(x_prior,"x_prior");
  g.add_node(t_mean,"t_mean");
  g.add_node(t_sigma,"t_sigma");
  g.add_node(t_dof,"t_dof");
  g.add_node(x_likelihood,"x_likelihood",{{t_mean,0},{t_sigma,0},{t_dof,0},{x_prior,0}});
  g.initialize();
  auto mx=g.get_monitor("x_prior",0);
  
  std::vector<double> buffer;
  for(int i=0;i<100000;++i)
    {
      g.sample(rnd);
      //cout<<mx()<<endl;
      buffer.push_back(mx());
    }
  std::sort(buffer.begin(),buffer.end());
  cout<<setprecision(10);
  for(int i=0;i<buffer.size();++i)
    {
      cout<<buffer[i]<<" "<<i<<endl;
    }
}
