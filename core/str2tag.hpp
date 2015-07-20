#ifndef STR2TAG
#define STR2TAG
#include <string>
#include <regex>
#include "tag.hpp"
#include "mcmc_exception.hpp"

namespace mcmc_utilities
{

  tag_t<std::string> tag(const char* cs)
  {
    std::string s(cs);
    std::regex pname_regex("([[a-zA-Z]+[a-zA-Z0-9_]*])");
    //std::regex idx_regex("([+-]?[0-9]+)");
    auto i=std::regex_token_iterator<std::string::iterator>(s.begin(),s.end(), pname_regex);
    std::regex_token_iterator<std::string::iterator> rend;
    if(i==rend)
      {
	throw mcmc_exception("no tag can be formed");
      }
    string pname=*i;
    std::vector<int> idx;
    i++;
    while(i!=rend)
      {
	idx.push_back(std::stoi(*i));
      }
    return tag_t<std::string>(pname,idx);
  }
}


#endif
