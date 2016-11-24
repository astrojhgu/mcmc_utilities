#ifndef PTRANSFORM_HPP
#define PTRANSFORM_HPP
/*
 *this file contains different implementations to the transform operations include thread based and openmp based
 */


#include <algorithm>
#include <vector>
#include <thread>


namespace mcmc_utilities
{
  template< class InputIt, class OutputIt, class UnaryOperation >
  OutputIt ptransform( InputIt first1, InputIt last1, OutputIt d_first,
		       UnaryOperation unary_op )
  {
    size_t nelem=std::distance(first1,last1);
    std::vector<std::thread> pool;
    for(auto i=first1;i!=last1;++i,++d_first)
      {
	auto x=*i;	
	pool.push_back(std::thread([&unary_op](const decltype(x) y,OutputIt it){*it=unary_op(y);},x,d_first));
      }
    std::for_each(pool.begin(),pool.end(),[&](std::thread& th){th.join();});
    return d_first;
  }
  
  template <typename T_vec1,typename T_vec2,typename UnaryOperation>
  void pvec_transform_thread(const T_vec1& vec_in,
			     T_vec2& vec_out,UnaryOperation unary_op,size_t max_th)
  {
    size_t nth=std::min(vec_in.size(),static_cast<decltype(vec_in.size())>(std::thread::hardware_concurrency()));
    if(nth>max_th)
      {
	nth=max_th;
      }
    size_t block_size=vec_in.size()/nth;    
    if(block_size*nth<vec_in.size())
      {
	block_size+=1;
      }
    std::vector<std::thread> pool;
    for(size_t i=0;i<nth;++i)
      {
	size_t nbegin=i*block_size;
	size_t nend=(i+1)*block_size;
	if(nend>=vec_in.size())
	  {
	    nend=vec_in.size();
	  }
	pool.push_back(std::thread([&vec_in,&vec_out,&unary_op](size_t nbegin,size_t nend){
	      for(size_t i=nbegin;i!=nend;++i)
		{
		  vec_out[i]=unary_op(vec_in[i]);
		}
	      for(;;);},nbegin,nend));
      }
    std::for_each(pool.begin(),pool.end(),[&](std::thread& th){if(th.joinable()) th.join();});
  }
  
  
  template <typename T_vec1,typename T_vec2,typename UnaryOperation>
  void pvec_transform_omp(const T_vec1& vec_in,
			  T_vec2& vec_out,UnaryOperation unary_op,size_t max_th)
  {
#pragma omp parallel for
    for(size_t i=0;i<vec_in.size();++i)
      {
	vec_out[i]=unary_op(vec_in[i]);
      }
  }
  
}


#endif
