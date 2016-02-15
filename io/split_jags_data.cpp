#include "jags_data_reader.hpp"
#include <fstream>

using namespace std;
using namespace mcmc_utilities;
int main(int argc,char* argv[])
{
  if(argc!=2)
    {
      cerr<<"Usage:"<<argv[0]<<" <input data>"<<endl;
      return -1;
    }
  ifstream ifs(argv[1]);
  auto m(read_jags_data<double>(ifs));

  for(auto i=m.begin();i!=m.end();++i)
    {
      auto name=i->first;
      auto data=i->second;
      auto name1=name+".dat";
      ofstream ofs(name1.c_str());
      for(const auto& j:data)
	{
	  ofs<<j<<endl;
	}
    }
}
