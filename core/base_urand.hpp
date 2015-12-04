#ifndef BASE_URAND
#define BASE_URAND

namespace mcmc_utilities
{
  template <typename T>
  class base_urand
  {
  public:
    T operator()()const
    {
      return do_rand();
    }
  private:
    virtual T do_rand()const=0;
      
  };
  

}

#endif
