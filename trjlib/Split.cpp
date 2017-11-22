/*
 * Split.cpp
 *
 *  Created on: Aug 8, 2017
 *      Author: marchi
 */

#include "Split.h"

#ifdef __INTEL
std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}
#else
std::vector<std::string> split(const std::string & s0, const std::string delim) {
	cout << "ukka"<<endl;

	vector<string> result;
	string s{s0};
	std::regex e{delim};

	std::regex_token_iterator<std::string::iterator> rend;
	std::regex_token_iterator<std::string::iterator> d ( s.begin(), s.end(), e, -1 );
	while (d!=rend) {
		string str=*d++;
		if(str.size() == 0) continue;
		result.push_back(str);
	}
	return result;
}
#endif
