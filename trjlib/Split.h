/*
 * Split.h
 *
 *  Created on: Aug 8, 2017
 *      Author: marchi
 */

#ifndef SPLIT_H_
#define SPLIT_H_

#include <iostream>
#include <string>
#include <vector>
#include <regex>

using std::string;
using std::vector;
vector<std::string> split(const std::string &, std::string);
string rmComments(const string & , const string );

#endif /* SPLIT_H_ */
