#include "mcmc.hpp"
#include <core/urand.hpp>
#include <cstdint>
namespace mcmc_utilities
{
  using plist_t=std::vector<std::pair<const node<double,std_vector>*,size_t> >;
  gtype* new_graph()
  {
    gtype* p=new gtype();
    std::cerr<<"new graph created @"<<(int64_t)p<<std::endl;
    return p;
  }

  void delete_graph(gtype* pg)
  {
    std::cerr<<"graph @"<<(int64_t)pg<<"destroied"<<std::endl;
    delete pg;
  }

  void sample(gtype* pg)
  {
    static urand<double> rnd;
    pg->sample(rnd);
  }

  void initialize(gtype* pg)
  {
    pg->initialize();
  }

  double get_value(const node<double,std_vector>* pn,int idx)
  {
    return pn->value(idx);
  }

  stochastic_node<double,std_vector>* set_observed_value(stochastic_node<double,std_vector>* pn,int idx,double v)
  {
    pn->set_observed(idx,true);
    pn->set_value(idx,v);
    return pn;
  }

  stochastic_node<double,std_vector>* set_value(stochastic_node<double,std_vector>* pn,int idx,double v)
  {
    pn->set_value(idx,v);
    return pn;
  }

  abs_node<double,std_vector>* add_abs_node(const node<double,std_vector>* x,int idx,gtype* pg)
  {
    auto pn=new abs_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,idx}});
    return pn;
  }

  add_node<double,std_vector>* add_add_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new add_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  sub_node<double,std_vector>* add_sub_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new sub_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  mul_node<double,std_vector>* add_mul_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new mul_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  div_node<double,std_vector>* add_div_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new div_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  pow_node<double,std_vector>* add_pow_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new pow_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  beta_node<double,std_vector>* add_beta_node(double a,double b,
					      gtype* pg)
  {
    auto pn=new beta_node<double,std_vector>(a,b);
    pg->add_node(pn,pn);
    return pn;
  }

  bin_node<double,std_vector>* add_bin_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new bin_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  bvnormal_node<double,std_vector>* add_bvnormal_node(const node<double,std_vector>* m1,int mi1,
						      const node<double,std_vector>* m2,int mi2,
						      const node<double,std_vector>* s1,int si1,
						      const node<double,std_vector>* s2,int si2,
						      const node<double,std_vector>* rho,int rhoi,
						      gtype* pg)
  {
    auto pn=new bvnormal_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{m1,mi1},{m2,mi2},{s1,si1},{s2,si2},{rho,rhoi}});
    return pn;
  }

  
  lt_node<double,std_vector>* add_lt_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new lt_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  le_node<double,std_vector>* add_le_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new le_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  gt_node<double,std_vector>* add_gt_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new gt_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  ge_node<double,std_vector>* add_ge_node(const node<double,std_vector>* x,int i1,
					   const node<double,std_vector>* y,int i2,
					   gtype* pg)
  {
    auto pn=new ge_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{x,i1},{y,i2}});
    return pn;
  }

  cond_node<double,std_vector>* add_cond_node(const node<double,std_vector>* s,int is,
					      const node<double,std_vector>* x,int ix,
					      const node<double,std_vector>* y,int iy,					      
					      gtype* pg)
  {
    auto pn=new cond_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{s,is},{x,ix},{y,iy}});
    return pn;
  }

  const_node<double,std_vector>* add_const_node(double v,gtype* pg)
  {
    auto pn=new const_node<double,std_vector>(v);
    pg->add_node(pn,pn);
    return pn;
  }

  cos_node<double,std_vector>* add_cos_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new cos_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  sin_node<double,std_vector>* add_sin_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new sin_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  tan_node<double,std_vector>* add_tan_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new tan_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  discrete_uniform_node<double,std_vector>* add_discrete_uniform_node(int a,int b,
					      gtype* pg)
  {
    auto pn=new discrete_uniform_node<double,std_vector>(a,b);
    pg->add_node(pn,pn);
    return pn;
  }

  gamma_node<double,std_vector>* add_gamma_node(const node<double,std_vector>* r,
					       int ri,
					       const node<double,std_vector>* l,
					       int li,
					       gtype* pg)
  {
    auto pn=new gamma_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{r,ri},{l,li}});
    return pn;
  }

  ilogit_node<double,std_vector>* add_ilogit_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new ilogit_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  logit_node<double,std_vector>* add_logit_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new logit_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  log_node<double,std_vector>* add_log_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new log_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  log10_node<double,std_vector>* add_log10_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new log10_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  normal_node<double,std_vector>* add_normal_node(const node<double,std_vector>* m,int midx,
						  const node<double,std_vector>* s,int sidx,
						   gtype* pg)
  {
    auto pn=new normal_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{m,midx},{s,sidx}});
    return pn;
  }

  pareto_node<double,std_vector>* add_pareto_node(const node<double,std_vector>* c,int cidx,
						  const node<double,std_vector>* a,int aidx,
						   gtype* pg)
  {
    auto pn=new pareto_node<double,std_vector>;
    pg->add_node(pn,pn,plist_t{{c,cidx},{a,aidx}});
    return pn;
  }

  phi_node<double,std_vector>* add_phi_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new phi_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  poisson_node<double,std_vector>* add_poisson_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new poisson_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  sqrt_node<double,std_vector>* add_sqrt_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new sqrt_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  step_node<double,std_vector>* add_step_node(const node<double,std_vector>* x,int i,gtype* pg)
  {
    auto pn=new step_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{x,i}});
    return pn;
  }

  t_node<double,std_vector>* add_t_node(const node<double,std_vector>* mu,int imu,
					const node<double,std_vector>* sigma,int isigma,
					const node<double,std_vector>* tau,int itau,
					gtype* pg)
  {
    auto pn=new t_node<double,std_vector>();
    pg->add_node(pn,pn,plist_t{{mu,imu},{sigma,isigma},{tau,itau}});
    return pn;
  }

  uniform_node<double,std_vector>* add_uniform_node(double a,double b,
					      gtype* pg)
  {
    auto pn=new uniform_node<double,std_vector>(a,b);
    pg->add_node(pn,pn);
    return pn;
  }

}

