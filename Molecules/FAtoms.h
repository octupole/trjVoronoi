/*
 * FAtoms.h
 *
 *  Created on: May 13, 2012
 *      Author: marchi
 */

#ifndef FATOMS_H_
#define FATOMS_H_
#include <iostream>
#include <string>
#include "Atoms.h"
#include "HeaderTrj.h"
#include "Ftypedefs.h"
#include "MyUtilClass.h"
//#include "nrjac.h"
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include "myEnums.hpp"

using namespace DVECT;
using namespace std;
typedef float Real;
using namespace Enums;
template <typename T>
class bambi: public Atoms<T> {

};
template <typename T>
class FAtoms: public Atoms<T> {
	using Atoms<T>::nr;
	using Atoms<T>::setCoord;
	using Atoms<T>::doCOtoOC;
	using Atoms<T>::x;
	using Atoms<T>::xa;
	using Atoms<T>::Mt;
	using Atoms<T>::rd;
protected:

public:
	using Atoms<T>::Atoms;
	virtual ~FAtoms();
	void setrd(vector<double> & rdx) {rd=rdx;};
	void setrd(Topol_NS::Topol & y){rd=y.getrd();};
	double getrd(int n){return rd[n];}
	void pdb(const vector<string> & );
};

#endif /* FATOMS_H_ */
