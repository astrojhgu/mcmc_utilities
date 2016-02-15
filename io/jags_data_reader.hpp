#ifndef JAGS_DATA_READER
#define JAGS_DATA_READER

#include <regex>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cassert>
namespace mcmc_utilities
{
  template <typename T>
  std::map<std::string,
	   std::vector<T> > read_jags_data(std::istream& is)
  {
    std::map<std::string,std::vector<T> > result;
    std::string input;
    std::string s;
    while(is>>s)
      {
	input+=s;
      }
    
    std::regex e1("\\s*<-\\s*c\\s*\\(");
    std::regex e2("\\s*\\)\\s*");
    std::regex e3("\\s*,\\s*");
    input=std::regex_replace(input,e1,std::string(" "));
    input=std::regex_replace(input,e2,std::string(" "));
    input=std::regex_replace(input,e3,std::string(" "));
    //std::cout<<input<<std::endl;
    std::istringstream iss(input);
    std::string name;
    std::vector<T> data;
    while(1)
      {
	std::string s;
	iss>>s;
	if(!iss.good())
	  {
	    break;
	  }
	std::istringstream iss1(s+'\n');
	T x;
	iss1>>x;
	//std::cout<<x<<std::endl;
	if(!iss1.good())
	  {
	    if(name!="")
	      {
		//std::cout<<name<<std::endl;
		result[name]=data;
	      }
	    name=s;
	    data.clear();
	  }
	else
	  {
	    data.push_back(x);
	  }
      }
    if(name!="")
      {
	//std::cout<<name<<std::endl;
	result[name]=data;
      }
    
    return result;
  }


  template <typename T>
  std::map<std::string,
	   std::vector<T> > read_jags_output(std::istream& is_index,
					   std::istream& is_data)
  {
    std::map<std::string,std::vector<T> > result;
    size_t current_line=0;
    for(;;)
      {
	std::string var;
	size_t line1,line2;

	is_index>>var>>line1>>line2;
	if(!is_index.good())
	  {
	    break;
	  }
	std::vector<T> data;
	for(;current_line<line2;++current_line)
	  {
	    std::string line;
	    std::getline(is_data,line);
	    if(!is_data.good())
	      {
		break;
	      }
	    line+=" ";
	    std::istringstream iss(line);
	    int dummy;
	    T d;
	    iss>>dummy>>d;
	    if(!iss.good())
	      {
		break;
	      }
	    data.push_back(d);
	  }
	result[var]=data;
      }
    return result;
  }
}

#endif
