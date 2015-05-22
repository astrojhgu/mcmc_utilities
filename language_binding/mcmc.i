%module(directors="1",allprotected="1") mcmc
%include "std_vector.i";
%include "stl.i"
%include "exception.i"

%{
#define private public
#include "core/gibbs_sampler.hpp"
#include "core/distribution.hpp"
#include "core/arms.hpp"
#include "core/mcmc_exception.hpp"
#include "vwrap.hpp"
#include "math/distributions.hpp"
#include "math/functions.hpp"
  using namespace std;
  using namespace mcmc_utilities;
  %}

%feature("autodoc", "1");
%feature("director") vector<double>;
%feature("director") dist_md<double,std::vector<double> >;
%apply double {std::vector<double,std::allocator<double > >::value_type };
%apply double {std::vector<double >::value_type };
%apply double {element_type_trait<std::vector<double> >::element_type};

%apply const double& {const std::vector<double,std::allocator<double > >::value_type& };
%apply const double& {const std::vector<double >::value_type& };

//%template (bpdensity) mcmc_utilities::probability_density_md<double,std::vector<double> >;
%include "vwrap.hpp"

#define private public
%include "core/gibbs_sampler.hpp"
%include "core/distribution.hpp"
%include "core/arms.hpp"
%include "core/mcmc_exception.hpp"
%include "math/distributions.hpp"
%include "math/functions.hpp"

%exception
{
  try
    {
      $action;
    }
  catch(const mcmc_exception& e)
    {
      std::cerr<<"mcmc error:"<<e.what()<<endl;
      SWIG_exception(SWIG_RuntimeError,e.what());
    }
}

%apply SWIGTYPE EXCEPTION_BY_VAL {mcmc_exception};

%template (array) std::vector<double>;
%template (gibbs_sample_real) mcmc_utilities::gibbs_sample<double,std::vector<double> , mcmc_utilities::u_random<double> >;
%template (pdensity) mcmc_utilities::dist_md<double,std::vector<double> >;
%template (urand) mcmc_utilities::u_random<double>;

%template (logdgamma) mcmc_utilities::logdgamma<double,double>;
%template (logdt) mcmc_utilities::logdt<double,double>;
%template (logdnorm) mcmc_utilities::logdnorm<double,double>;
%template (logdlnorm) mcmc_utilities::logdlnorm<double,double>;
%template (logdbin) mcmc_utilities::logdbin<double,double>;
%template (logdpoisson) mcmc_utilities::logdpoisson<double,double>;
%template (logdpar) mcmc_utilities::logdpar<double,double>;
