#ifndef BASE_URAND
#define BASE_URAND

namespace mcmc_utilities
{
  template <typename T>
  class base_urand
  {
  public:
    virtual ~base_urand(){}
  public:
    T operator()()
    {
      return do_rand();
    }

  public:
    T max()const
    {
      return T(1);
    }

    T min()const
    {
      return T(0);
    }

    bool is_parallel()const
    {
      return do_is_parallel();
    }
  private:
    virtual T do_rand()=0;
    virtual bool do_is_parallel()const
    {
      return false;
    }
      
  };
  

}

#endif
