#include "jags_data_reader.hpp"
#include <fstream>

using namespace std;
using namespace mcmc_utilities;
int main(int argc,char* argv[])
{
  if(argc!=3)
    {
      cerr<<"Usage:"<<argv[0]<<" <input index file> <input data>"<<endl;
      return -1;
    }
  ifstream ifs_index(argv[1]);
  ifstream ifs_data(argv[2]);
  auto m(read_jags_output<double>(ifs_index,ifs_data));

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
