#include <core/graph.hpp>
#include <math/distributions.hpp>
#include <core/stochastic_node.hpp>
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
	result+=logdbin(value(0,i),parent(0,i),parent(1,i));

      }
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

  g.add_node(new uniform(0,1),"A",{});//vnode("A")=uniform(0,1)
  g.add_node(new uniform(0,1),"B",{});//vnode("B")=uniform(0,1)

  g.add_node(new uniform(0,100),"mu",{});//vnode("mu")=uniform(0,100)
  g.add_node(new uniform(0,100),"sigma",{});//vnode("sigma")=uniform(0,100)

  g.add_node(new eff(),"eff",{{"A",0},{"B",0},{"E",0},{"mu",0},{"sigma",0}});
  //vnode("eff")=vnode("A")+(vnode("B")-vnode("A"))*phi((vnode("E")-vnode("mu"))/vnode("sigma"));
  

  g.add_node(dl.get_nrec(),"nrec",{{"eff",0},{"ninj",0}});
  //vnode("nrec")=bin(vnode("eff",0),vnode("ninj",0))
  
  g.set_value({"mu"},0,30.0);

  for(int i=0;i<100;++i)
    {
      g.sample(rnd1);
    }
  
  for(int i=0;i<30000;++i)
    {
      g.sample(rnd1);
      auto p(g.get_params());
      for(auto& x:p)
	{
	  cout<<x<<" ";
	}
      cout<<endl;
    }
}
