/*
 * VoronoiBinary.cpp
 *
 *  Created on: Dec 2, 2017
 *      Author: marchi
 */

#include "VoronoiBinary.h"

namespace Voro {
void VoronoiBinary::WriteIt(std::ofstream & fout){
	if(writeBinary)
		bPrintBody(fout);
	else
		VoronoiMicelles::WriteIt(fout);

}
void VoronoiBinary::ReadIt(std::ifstream & fin){
	bReadBody(fin);
}

VoronoiBinary::VoronoiBinary(ifstream & fin) {
	bReadHeader(fin);
	Vol=vector<double>(nr);
	Neighs=vector<vector<int>>(nr);
	Surface=vector<vector<double>>(nr);
	area.Allocate(nresid,nc);
	Vols.Allocate(nresid);
	wShells=vector<vector<int>>(VoronoiSetter::maxLevel);

	this->readBinary=true;
}
VoronoiBinary::VoronoiBinary(ofstream & fout,Topol & myTop,bool bH): VoronoiMicelles::VoronoiMicelles(myTop,bH){
	this->bPrintHeader(fout);
	this->writeBinary=true;
}

VoronoiBinary::~VoronoiBinary() {
	// TODO Auto-generated destructor stub
}

} /* namespace Voro */
