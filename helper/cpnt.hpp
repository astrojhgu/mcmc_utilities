#ifndef CPNT_HPP
#define CPNT_HPP
#include <memory>


namespace mcmc_utilities
{
  template <typename T_p,typename T_var1>
  class cpnt
  {
  public:
    std::shared_ptr<_vnode<T_p,T_var1> > pn;
    size_t n;
    
  public:
    cpnt(const std::shared_ptr<_vnode<T_p,T_var1> > _pn,
	 size_t _n)
      :pn(_pn),n(_n)
    {}
    
    cpnt(const std::shared_ptr<_vnode<T_p,T_var1> > _pn)
      :cpnt(_pn,0)
    {}

  public:
    operator std::shared_ptr<_vnode<T_p,T_var1> >()
    {
      return pn;
    }
  };
}
#endif
