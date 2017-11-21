/*
 * Split.cpp
 *
 *  Created on: Aug 8, 2017
 *      Author: marchi
 */

#include "Split.h"

std::vector<std::string> split(const std::string & s0, const std::string delim) {
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
string rmComments(const string & s, const string comment){
	std::regex e{comment};
	string result{std::regex_replace(s,e,"")};
	return result;
}
