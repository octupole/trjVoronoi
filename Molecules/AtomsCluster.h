/*
 * AtomsCluster.h
 *
 *  Created on: Nov 28, 2017
 *      Author: marchi
 */

#ifndef MOLECULES_ATOMSCLUSTER_H_
#define MOLECULES_ATOMSCLUSTER_H_

#include "FAtoms.h"
#include "Percolation.h"
#include "TopolMic.h"
#include <algorithm>

using namespace Topol_NS;
using std::vector;
using std::min_element;

template <typename T>
class AtomsCluster: public FAtoms<T> {
	static int calls;
	const T HALF{0.5000-0.0001};
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	Percolation<T> * Perco{nullptr};
	void __ReconstructOneCluster(vector<bool> &);
	Dvect __FindCell(const vector<vector<int>> & ,const vector<vector<int>> & );
	using Atoms<T>::x;
	using Atoms<T>::xa;
	using Atoms<T>::Mt;
	using Atoms<T>::rd;
	using Atoms<T>::nr;
public:
	using FAtoms<T>::FAtoms;
	void Reconstruct(const string &, TopolMic &);
	void Reconstruct(Contacts<T> *);
	void SetupPercolate();
	void Percolate();
	void Percolate(double y){
		Perco->setRcut(y);
		this->Percolate();
	}
	virtual bool hasPerco() const {if(Perco) return true;return false;}

	virtual ~AtomsCluster();
};

#endif /* MOLECULES_ATOMSCLUSTER_H_ */
