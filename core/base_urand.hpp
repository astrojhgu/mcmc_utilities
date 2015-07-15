#ifndef BASE_URAND
#define BASE_URAND

namespace mcmc_utilities
{
  template <typename T_p>
  class base_urand
  {
  public:
    T_p operator()()const
    {
      return do_rand();
    }
  private:
    virtual T_p do_rand()const=0;
      
  };
  

}

#endif
