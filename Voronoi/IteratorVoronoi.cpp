/*
 * IteratorMAtoms.cpp
 *
 *  Created on: Dec 17, 2015
 *      Author: marchi
 */

#include "IteratorVoronoi.h"

namespace myiterators {

IteratorVoronoi::IteratorVoronoi(pointer x, ifstream * y,long int nstrt,long int nnd): p{x},
finx{y},nstart{nstrt}, nend{nnd}{
	cout << nstart << " " <<nend<<endl;
	cout <<p->getCO()<<endl;
	ios::streampos whereIwas=finx->tellg();
	finx->seekg(0,finx->end);
	cout << " Illa "<< whereIwas << " " <<endl;
	len=finx->tellg();
	finx->seekg(0,finx->beg);
	finx->seekg(whereIwas,finx->beg);
	cout << " Illa "<< finx->tellg() << " " <<endl;
	exit(1);
//	len=finx->tellg()-whereIwas+1;


}


IteratorVoronoi & IteratorVoronoi::operator++(){
	ifstream & fin=*finx;
	while(finx->tellg() < len){
		if(ntime < nstart) {fin>> *p;ntime++;;continue;}
		if(nend != -1 && ntime > nend) break;

		fin >> *p;

		ntime++;
		return *this;
	}
	(*this)=IteratorVoronoi();
	return *this;
}



} /* namespace myiterators */
