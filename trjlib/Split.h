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
#ifdef __INTEL
vector<std::string> split(const std::string &, char);
#else
vector<std::string> split(const std::string &, std::string);
#endif
#endif /* SPLIT_H_ */
