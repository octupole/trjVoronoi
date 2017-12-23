/*
 * PercolationJSON.h
 *
 *  Created on: Dec 21, 2017
 *      Author: marchi
 */

#ifndef MOLECULES_PERCOLATIONJSON_H_
#define MOLECULES_PERCOLATIONJSON_H_

#include "Percolation.h"
#include "json.hpp"
using nlohmann::json;

template <typename T>
class PercolationJSON: public Percolation<T> {
	json myJson;
	void __Writeit(ostream &);
public:
	using Percolation<T>::Percolation;
	virtual ~PercolationJSON();
};

#endif /* MOLECULES_PERCOLATIONJSON_H_ */
