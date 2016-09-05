#ifndef MCMC_HPP
#define MCMC_HPP
#include <vector>
#include <core/tag_t.hpp>
#include <core/graph.hpp>
#include <nodes/abs_node.hpp>
#include <nodes/arithmetic_node.hpp>
#include <nodes/beta_node.hpp>
#include <nodes/bin_node.hpp>
#include <nodes/bvnormal_node.hpp>
#include <nodes/compare_node.hpp>
#include <nodes/cond_node.hpp>
#include <nodes/const_node.hpp>
#include <nodes/cos_node.hpp>
#include <nodes/discrete_uniform_node.hpp>
#include <nodes/gamma_node.hpp>
#include <nodes/ilogit_node.hpp>
#include <nodes/log10_node.hpp>
#include <nodes/logit_node.hpp>
#include <nodes/log_node.hpp>
#include <nodes/normal_node.hpp>
#include <nodes/pareto_node.hpp>
#include <nodes/phi_node.hpp>
#include <nodes/poisson_node.hpp>
#include <nodes/sin_node.hpp>
#include <nodes/sqrt_node.hpp>
#include <nodes/step_node.hpp>
#include <nodes/switch_node.hpp>
#include <nodes/tan_node.hpp>
#include <nodes/t_node.hpp>
#include <nodes/trunc_pareto_node.hpp>
#include <nodes/uniform_node.hpp>



namespace mcmc_utilities
{
  template <typename T>
  using std_vector=std::vector<T>;

  using gtype=graph<double,const node<double,std_vector>*,std_vector>;
    
  gtype* new_graph();

  void delete_graph(gtype* pg);

  stochastic_node<double,std_vector>* set_observed_value(stochastic_node<double,std_vector>* pn,int idx,double v);
  
  stochastic_node<double,std_vector>* set_value(stochastic_node<double,std_vector>* pn,int idx,double v);  

  const_node<double,std_vector>* add_const_node(double v,gtype* pg);

  normal_node<double,std_vector>* add_normal_node(const node<double,std_vector>* m,int midx,
						  const node<double,std_vector>* s,int sidx,
						  gtype* pg);

  void sample(gtype* p);
}

#endif
