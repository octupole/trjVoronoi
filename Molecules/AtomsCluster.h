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
class CellComp{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	Dvect Xd;
	Matrix co;
public:
	CellComp(Dvect& x, Dvect & y, Matrix & c): Xd{x-y}, co{c}{};
	CellComp(Dvect& x, Matrix & c): Xd(x), co(c){};
	bool operator()(const Dvect & x, const Dvect & y){
		Dvect x1=Xd+x;
		Dvect x2=Xd+y;
		Dvect xc1=co*x1;
		Dvect xc2=co*x2;
		return xc1.Norm() < xc2.Norm();
	}
};
template <typename T>
class PBCvect{
	using Dvect=DDvect<T>;
	static vector<Dvect> * nxyz;
	const int M;
public:
	PBCvect(): M(3){
		nxyz=new vector<Dvect>;
		for(int o=0;o<M;o++){
			double nx=static_cast<double>(o-1);
			for(int p=0;p<M;p++){
				double ny=static_cast<double>(p-1);
				for(int q=0;q<M;q++){
					double nz=static_cast<double>(q-1);
					nxyz->push_back(Dvect(nx,ny,nz));
				}
			}
		}
	}
	static vector<Dvect> & getVec(){return *nxyz;}
};


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
	virtual void SetupPercolate(Topol_NS::Topol &);
	void SetupPercolate();
	int Percolate();
	int Percolate(double y){
		Perco->setRcut(y);
		return this->Percolate();
	}
	Percolation<T> * gPerco(){return Perco;}
	virtual bool hasPerco() const {if(Perco) return true;return false;}

	virtual ~AtomsCluster();
};

#endif /* MOLECULES_ATOMSCLUSTER_H_ */
