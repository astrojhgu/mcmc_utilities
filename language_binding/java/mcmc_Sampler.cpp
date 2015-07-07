#include "mcmc_Sampler.h"
#include <iostream>
#include <core/gibbs_sampler.hpp>
#include <vector>
#include <cassert>
using namespace std;
using namespace mcmc_utilities;

u_random<double> rng;

class JSampleable
  :public probability_density_md<double,std::vector<double> >
{
public:
  JNIEnv* env;
  jobject* p_ISampleable;
  jmethodID eval_id;
  jmethodID range_id;
  jmethodID ip_id;
  jmethodID set_id;
  jmethodID get_id;
  jmethodID ssize_id;
  jmethodID gsize_id;
  jmethodID pinit_id;
  jclass param_class;
public:
  JSampleable(JNIEnv* env1,
	      jobject* p_is,
	      jclass pclass)
    :env(env1),p_ISampleable(p_is),eval_id(0),range_id(0),ip_id(0),set_id(0),get_id(0),ssize_id(0),gsize_id(0),pinit_id(0),param_class(pclass)
  {
    jclass clazz=env->GetObjectClass(*p_ISampleable);
    eval_id=env->GetMethodID(clazz,"evalLog","(Lmcmc/IParameter;)D");
    range_id=env->GetMethodID(clazz,"varRange","(Lmcmc/IParameter;I)[D");
    ip_id=env->GetMethodID(clazz,"initPoints","(Lmcmc/IParameter;I)[D");
    get_id=env->GetMethodID(param_class,"getValue","(I)D");
    set_id=env->GetMethodID(param_class,"setValue","(ID)V");
    gsize_id=env->GetMethodID(param_class,"getSize","()I");
    ssize_id=env->GetMethodID(param_class,"setSize","(I)V");
    pinit_id=env->GetMethodID(param_class,"<init>","()V");
    assert(eval_id!=0);
    assert(range_id!=0);
    assert(ip_id!=0);
    assert(set_id!=0);
    assert(get_id!=0);
    assert(gsize_id!=0);
    assert(ssize_id!=0);
    assert(pinit_id!=0);
  }
public:
  jobject create_jparam(const std::vector<double>& x)const
  {
    jobject result=env->NewObject(param_class,pinit_id);
    env->CallVoidMethod(result,ssize_id,x.size());

    for(int i=0;i<x.size();++i)
      {
	env->CallVoidMethod(result,set_id,i,x[i]);
      }
    return result;
  }

  std::vector<double> param2vec(jobject o)
  {
    int length=env->CallIntMethod(o,gsize_id);
    std::vector<double> result(length);
    for(int i=0;i<length;++i)
      {
	result[i]=env->CallDoubleMethod(o,get_id,i);
      }
    return result;
  }
private:
  
  double do_eval_log(const std::vector<double>& x)const
  {
    jobject p=create_jparam(x);
    jdouble result=env->CallDoubleMethod(*p_ISampleable,eval_id,p);
    env->DeleteLocalRef(p);
    return result;
  }

  std::pair<double,double> do_var_range(const std::vector<double>& x,size_t ndim)const
  {
    jobject p=create_jparam(x);
    jdoubleArray result=static_cast<jdoubleArray>(env->CallObjectMethod(*p_ISampleable,range_id,p,ndim));
    jboolean isCopy=false;
    jdouble* data=env->GetDoubleArrayElements(result,&isCopy);
    jsize length=env->GetArrayLength(result);

    return make_pair(data[0],data[1]);
    env->ReleaseDoubleArrayElements(result,data,0);
    env->DeleteLocalRef(p);
  }

  std::vector<double> do_init_points(const std::vector<double>& x,size_t ndim)const
  {
    jobject p=create_jparam(x);
    jdoubleArray result=static_cast<jdoubleArray>(env->CallObjectMethod(*p_ISampleable,ip_id,x,ndim));
    jboolean isCopy=false;
    jdouble* data=env->GetDoubleArrayElements(result,&isCopy);
    jsize length=env->GetArrayLength(result);

    //return make_pair(data[0],data[1]);
    std::vector<double> result1(length);

    for(int i=0;i<length;++i)
      {
	result1[i]=data[i];
      }
    env->ReleaseDoubleArrayElements(result,data,0);
    env->DeleteLocalRef(p);
    return result1;
  }
};



JNIEXPORT void JNICALL Java_mcmc_Sampler_gibbsSample
(JNIEnv * env, jclass , jobject is, jobject init_var, jboolean do_metropolis)
{
  jclass clazz=env->GetObjectClass(init_var);
  
  JSampleable js(env,&is,clazz);

  std::vector<double> p(js.param2vec(init_var));
  
  gibbs_sample(js,p,1,rng);


  for(int i=0;i<p.size();++i)
    {
      env->CallVoidMethod(init_var,js.set_id,i,p[i]);
    }
}
