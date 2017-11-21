/*
 * FAtoms.cpp
 *
 *  Created on: May 13, 2012
 *      Author: marchi
 */

#include "FAtoms.h"


template <typename T>
void FAtoms<T>::pdb(const vector<string> & data){

	vector<Dvect> xc;
	double aa=-1.0,bb=-1.0,cc=-1.0,alpha=90.0,beta=90.0,gamma=90.0;
	for(size_t i=0;i<data.size();i++){
		if(data[i].find("CRYST1") == 0){

			istringstream ss(data[i].substr(6,49));
			ss>>aa>>bb>>cc>>alpha>>beta>>gamma;
			aa*=unit_nm;
			bb*=unit_nm;
			cc*=unit_nm;
		}
	}
	try{
	if(aa < 0 || bb < 0 || cc <0) throw string("\nWarning: PDB file does not have CRYST1 keyword. ")+
			string("If you are doing calculations\n     on this structure, box is not set.\n");
	} catch(const string & s){cout << s<<endl;};
	Metric<T> Met(brot<T>()(aa,bb,cc,alpha,beta,gamma));
	for(size_t i=0;i<data.size();i++){
		if(data[i].find("ATOM") == 0 || data[i].find("HETATM") == 0 ) {
			Dvect tmp;
			std::stringstream ss(data[i].substr(30,24));
			ss>> tmp;
			xc.push_back(tmp);
		}
	}
	try{
		if(nr != int(xc.size())) throw string("Something wrong: No. atoms do not match");
	} catch(const string & s){
		cout << s <<endl;
		exit(1);
	}

	rvec * x0=new rvec[nr];
	for(int n=0;n<nr;n++){
		for(int o=0;o<DIM;o++) x0[n][o]=xc[n][o]*unit_nm;
	}
	setCoord(Met,x0);
	doCOtoOC();
	delete [] x0;
}


template <typename T>
FAtoms<T>::~FAtoms() {}

template class FAtoms<float>;
template class FAtoms<double>;

