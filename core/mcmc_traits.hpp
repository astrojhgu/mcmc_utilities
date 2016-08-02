#ifndef MCMC_TRAITS
#define MCMC_TRAITS
#define MCMC_HEADER
#include <cmath>
#include <cstddef>
namespace mcmc_utilities
{
  /**
     Get the size of an array object
     \param x the array object
     \return the size of the array object
   */
  template <typename T>
  inline size_t get_size(const T& x)
  {
    return x.size();
  }
  
  /**
     \brief Trait class, in which the types of elements in an array are defined
     \tparam the type of the array object
   */
  template <typename T>
  class element_type_trait
  {
  public:
    /**
       Default definition of element_type
     */
    typedef typename T::value_type element_type;
  };

  
  /**
     \brief The return type trait of some certain data types.
  */
  template <typename T>
  class return_type_trait
  {
  public:
    typedef T value_type;
    typedef T& reference_type;
    typedef const T& const_reference_type;
  };
  

  /**
     Help function to get the i-th element from an array
     \tparam T the type of the array object
     \param x the array object
     \param i the order of the element
     \return the fetched element value, const reference
   */
  template <typename T>
  inline typename 
  return_type_trait<typename element_type_trait<T>::element_type>::
  const_reference_type get_element(const T& x,size_t i)
  {
    return x[i];
  }

  
  /**
     set ths i-th element by a given value
     \tparam T the type of the array object
     \tparam Tx the type of the element
     \param x the array object
     \param i the order of the element
     \param v the value of the element to be set
   */
  template<typename T,typename TX>
  inline void set_element(T& x,size_t i,
			  const TX& v)
  {
    x[i]=v;
  }

  template <typename T,typename TX>
  inline void push_back(T& x,const TX& v)
  {
    x.push_back(v);
  }

  template <typename T>
  inline typename return_type_trait<typename element_type_trait<T>::element_type>::
  const_reference_type last_element(const T& x)
  {
    return x.back();
  }

  template <typename T>
  inline typename return_type_trait<typename element_type_trait<T>::element_type>::
  reference_type last_element(T& x)
  {
    return x.back();
  }


  /**
     resize an array object
     \tparam T the type of the array
     \param x the array object
     \param s the new size
   */
  template <typename T>
  inline void resize(T& x,size_t s)
  {
    x.resize(s);
  }

  template <typename T>
  inline void reserve(T& x,size_t s)
  {
    x.reserve(s);
  }


  /**
     Assignment operator of two array objects
     \tparam Tl the type of left-hand array
     \tparam Tr the type of right-hand array
     \param lhs the left-hand array
     \param rhs the right-hand array
     \return the reference of the left-hand array
   */
  template <typename Tl,typename Tr>
  inline Tl& mcmc_assign(Tl& lhs,const Tr& rhs)
  {
    return (lhs=rhs);
  }

  template <typename T>
  constexpr T C_ONE()
  {
    return static_cast<T>(1);
  }

  template <typename T>
  constexpr T C_ZERO()
  {
    return static_cast<T>(0);
  }

  template <typename T,unsigned int N>
  constexpr T C_N2T()
  {
    return static_cast<T>(N);
  }
  
  template <typename T>
  constexpr T C_NAN()
  {
    return static_cast<T>(std::nan(""));
  }
}



#endif
