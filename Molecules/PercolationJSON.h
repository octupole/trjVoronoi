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

template <typename T>
class PercolationJSON: public Percolation<T> {
public:
	using Percolation<T>::Percolation;
	virtual ~PercolationJSON();
};

#endif /* MOLECULES_PERCOLATIONJSON_H_ */
