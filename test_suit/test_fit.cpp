#include <iostream>
#include <vector>
#include <fstream>
#include <fio.h>
#include "core/gibbs_sampler.hpp"
using namespace std;
using namespace mcmc_utilities;
using namespace std;
std::vector<double> x_vec;
std::vector<double> y_vec;

class std_norm
  :public probability_density_md<double,std::vector<double> >
{
  double do_eval_log(const std::vector<double>& param)const
  {
    double p=0;
    double k=param[0];
    double b=param[1];
    double sigma=param[2];
    
    for(int i=0;i<x_vec.size();++i)
      {
	double y1=k*x_vec[i]+b;
	p+=((-(y_vec[i]-y1)*(y_vec[i]-y1)/(2*sigma*sigma))-std::log(sigma));
      }
    
    //p=-k*k-b*b-sigma*sigma/4;
    //cout<<"a ";
    for(int i=0;i<3;++i)
      {
	//cout<<param[i]<<" ";
      }
    //cout<<p<<endl;
    
    return p;
  }

  probability_density_md* do_clone()const
  {
    return new std_norm(*this);
  }

  void do_var_range(std::vector<double>& x1,std::vector<double>& x2)const
  {
    for(int i=0;i<x1.size();++i)
      {
	x1[i]=-20;
	x2[i]=20;
      }
    x1[2]=.001;
  }

  
};


int main()
{
  ifstream ifs("data.txt");
  for(;;)
    {
      double x,y;
      ifs>>x>>y;
      if(!ifs.good())
	{
	  break;
	}
      x_vec.push_back(x);
      y_vec.push_back(y);
    }
  //uniform_rng<double,double> ur;
  std_norm sn;
  std::vector<double> init_var(3,0);
  init_var[0]=4;
  init_var[1]=3;
  init_var[2]=20;
  
  //cout<<sn.eval(init_var)<<endl;
  
  int s=(time(0));
  
  //cout<<"seed:"<<s<<endl;
  srand(s);
  //return 0;
  init_var[0]=4;
  init_var[1]=3;
  init_var[2]=.5;
  
  for(int i=0;i<100000;++i)
    {
      gibbs_sample(sn,init_var,1,u_random<double>);

      cout<<i<<" ";
      //cout<<init_var[1]<<" "<<init_var[2]<<endl;
      for(int j=0;j<init_var.size();++j)
	{
	  cout<<init_var[j]<<" ";
	}
      cout<<endl;
      
    }
}
