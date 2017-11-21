/*
 * LCells.h
 *
 *  Created on: Dec 18, 2014
 *      Author: marchi
 */

#ifndef LINKEDCELLS_SRC_LCELLS_H_
#define LINKEDCELLS_SRC_LCELLS_H_

#include <vector>
#include <algorithm>
#include "Ftypedefs.h"
#include "MyUtilClass.h"



using namespace std;
using namespace MATRIX;
using namespace DVECT;


template <typename T>
class LCells {
	const T HALF{0.5000-0.0001};
	const T small{1.0e-6};
	using Matrix=MMatrix<T>;
	using Dvect=DDvect<T>;
	typedef vector<int> vectint;
	int nr{0};
	T Rcut{0.0};
	const T Rmax{45};

	vector<int> nc={-1,-1,-1};
	vector<vectint>  indx;
	vector<vectint> Cellp;
	vector<vector<vector<vectint> > > Chainp;
	Matrix co,oc;
	vector<Dvect> x;
	vector<vector<int> > nnl;
	T Dist_ijk(int, int,int);
	void Init(Matrix & co0, const vector<Dvect> & y){
		T Rcut0=Rcut;
		x=y;
		co=co0;
		oc=co.Inversion();
		nr=x.size();

		try{
			if(co[XX][XX] < Rcut*2.0) throw string{"Box is too small to run with neighbor lists."};
		}catch(const string & s){
			cout << s<<endl;
			exit(1);
		}
		if(co[XX][XX] > Rmax){
			if(nc[XX] < 0){
				nc[XX]=static_cast<int>(co[XX][XX]/Rcut0);
				nc[YY]=static_cast<int>(co[YY][YY]/Rcut0);
				nc[ZZ]=static_cast<int>(co[ZZ][ZZ]/Rcut0);
			}
			if(Chainp.size()) Chainp.clear();
			Chainp=vector<vector<vector<vectint> > >(nc[XX]);
			for(int m=0;m<nc[XX];m++){
				Chainp[m]=vector<vector<vectint> >(nc[YY]);
				for(int n=0;n<nc[YY];n++)
					Chainp[m][n]=vector<vectint>(nc[ZZ]);
			}
		}
		if(nnl.size()) nnl.clear();
		nnl=vector<vectint>(nr);
	}
public:
	LCells();
	LCells(Matrix & co0, const vector<Dvect> & y, T rcut): Rcut{rcut}{
		this->Init(co0,y);
		this->Index();
	}
	void operator()(Matrix & co0, const vector<Dvect> y){
		this->Init(co0,y);
	}
	vector<vector<int> > & getNeigh(){
		try{
			if(!nnl.size()) throw " Neighbour list not yet computed !";
		}
		catch(const char * s){
			cout << s << endl;
			exit(1);
		}
		return nnl;
	}
	void Index();
	vector<vector<int> > &List(bool=true);
	virtual ~LCells();
};

#endif /* LINKEDCELLS_SRC_LCELLS_H_ */
