#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
#include <core/urand.hpp>
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

  double do_value(size_t idx,size_t obsid)const override
  {
    double A=this->parent(0,obsid);
    double B=this->parent(1,obsid);
    double E=this->parent(2,obsid);
    double mu=this->parent(3,obsid);
    double sigma=this->parent(4,obsid);
    return A+(B-A)*phi((E-mu)/sigma);
    //return this->parents[0]
  }
};


class uniform
  :public stochastic_node<double,double>
{
private:
  double a,b;
public:
  uniform(double a_,double b_)
    :stochastic_node<double,double>(0,(a_+b_)/2),a(a_),b(b_)
  {
  };

  double do_log_prior_prob()const override
  {
    return 0;
  }

  std::pair<double,double> do_var_range()const override
  {
    return std::make_pair(a,b);
  }

  bool is_continuous(size_t)const override
  {
    return true;
  }
};




class nrec
  :public observed_node<double,double>
{
public:
  std::vector<double> data;
public:
  nrec()
    :observed_node<double,double>(2)
  {}

  double do_log_prior_prob()const override
  {
    double result=0;
    for(int i=0;i<nobs();++i)
      {
	//result+=logdbin(value(0,i),parent(0,i),parent(1,i));
	double result1=logdbin(value(0,i),parent(0,i),parent(1,i));
	if(std::isnan(result1))
	  {
	    std::cerr<<value(0,i)<<" "<<parent(0,i)<<" "<<parent(1,i)<<endl;
	  }
	   
	result+=result1;
      }
    assert(!std::isnan(result));
    return result;
  }

  double do_value(size_t idx,size_t i)const override
  {
    return data.at(i);
  }

  size_t do_nobs()const override
  {
    return data.size();
  }


};

class ninj
  :public observed_node<double,double>
{
 public:
  std::vector<double> data;

  ninj()
    :observed_node<double,double>(0)
  {}

  double do_log_prior_prob()const override
  {
    assert(0);
    return 0;
  }

  double do_value(size_t idx,size_t i)const override
  {
    return data[i];
  }
  
  size_t do_nobs()const override
  {
    return data.size();
  }

};

class energy
  :public observed_node<double,double>
{
public:
  std::vector<double> data;
  energy()
    :observed_node<double,double>(0)
  {}

  double do_log_prior_prob()const override
  {
    assert(0);
    return 0;
  }

  double do_value(size_t idx,size_t i)const override
  {
    return data[i];
  }

  size_t do_nobs()const override
  {
    return data.size();
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

  shared_ptr<energy> get_energy()const
  {
    std::shared_ptr<energy> p(new energy());
    p->data=vec_e;
    return p;
  }

  shared_ptr<ninj> get_ninj()const
  {
    std::shared_ptr<ninj> p(new ninj());
    p->data=vec_ninj;
    return p;
  }

  shared_ptr<nrec> get_nrec()const
  {
    std::shared_ptr<nrec> p(new nrec());
    p->data=vec_nrec;
    return p;
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
  g.add_node(dl.get_energy(),"E",{});//vnode("E")=observed_constant()
  g.add_node(dl.get_ninj(),"ninj",{});//vnode("nij")=observed_constant()
  //g.add_node(dl.get_ninj(),{"nrec"},{});
  auto p=new uniform(0,1-1e-5);
  g.add_node(p,"A",{});//vnode("A")=uniform(0,1)
  g.add_node(new uniform(0,1-1e-5),"B",{});//vnode("B")=uniform(0,1)
  auto pmu=new uniform(0,100);
  g.add_node(pmu,"mu",{});//vnode("mu")=uniform(0,100)
  g.add_node(new uniform(0,100),"sigma",{});//vnode("sigma")=uniform(0,100)

  g.add_node(new eff(),"eff",{{"A",0},{"B",0},{"E",0},{"mu",0},{"sigma",0}});
  //vnode("eff")=vnode("A")+(vnode("B")-vnode("A"))*phi((vnode("E")-vnode("mu"))/vnode("sigma"));
  

  g.add_node(dl.get_nrec(),"nrec",{{"eff",0},{"ninj",0}});
  //vnode("nrec")=bin(vnode("eff",0),vnode("ninj",0))
  
  g.set_value("A",0,5.14989e-05);
  g.set_value("B",0,0.999406);
  g.set_value("mu",0,13.2949);
  //g.set_value("mu",0,31.7934);
  //g.set_value("sigma",0,17.5035);
  g.set_value("sigma",0,17.6002);

  //g.set_value("A",0,0.535597);
  //g.set_value("B",0,0.996047);

  
  auto A=g.get_monitor("A",0);
  auto B=g.get_monitor("B",0);
  auto mu=g.get_monitor("mu",0);
  auto sigma=g.get_monitor("sigma",0);
  double x=mu();
  for(int i=0;i<10000;++i)
    {
      g.sample(rnd1);
      //cout<<arms(*pmu,x,1,rnd1)<<endl;
      //if(i%100==0)
	{
	  cout<<A()<<" "<<B()<<" "<<mu()<<" "<<sigma()<<endl;
	}
    }
}
