#ifndef PRNG_HPP
#define PRNG_HPP
#include <atomic>
#include <mutex>
#include "prng_engine.hpp"

namespace mcmc_utilities
{
  template <typename T>
  class prng
    :public base_urand<T>
  {
  private:
    std::atomic<int> nused;
    std::mutex mtx;
    std::vector<sitmo::prng_engine> eng_array;
  public:
    prng()
      :eng_array(1),nused(0)
    {}

    prng(int n)
      :eng_array(n),nused(0)
    {
      for(int i=0;i<eng_array.size();++i)
	{
	  eng_array[i].seed(i);
	}
    }

  private:
    T do_rand()override
    {
      int n=nused++;
      if(n>=eng_array.size())
	{
	  mtx.lock();
	  eng_array.push_back(sitmo::prng_engine(n));
	  mtx.unlock();
	}
      //std::cerr<<n<<std::endl;
      T result=const_cast<sitmo::prng_engine&>(eng_array[n])()/static_cast<T>(sitmo::prng_engine::max());
      --nused;
      return result;
    }

    bool do_is_parallel()const override
    {
      return true;
    }
  };
}


#endif
