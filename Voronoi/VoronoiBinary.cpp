/*
 * VoronoiBinary.cpp
 *
 *  Created on: Dec 2, 2017
 *      Author: marchi
 */

#include "VoronoiBinary.h"

namespace Voro {
void VoronoiBinary::bReadBody(ifstream & fin){
	bool have_clusters{false};
	fin.read(as_byte(have_clusters),sizeof(have_clusters));
	fin.read(as_byte(nframe),sizeof(nframe));
	if(have_clusters){
		rdVector(fin,atClusters);
		rdVector(fin,VolClusters);
		rdVector(fin,Clusters);
		rdVector(fin,AreaClusters);
	}
	vector<double> Vol_;
	vector<vector<int>> Neighs_;
	vector<vector<double>> Surface_;
	rdVector(fin,Vol_);
	rdVector(fin,Neighs_);
	rdVector(fin,Surface_);
	for(size_t o{0};o<cindex.size();o++){
		Vol[cindex[o]]=Vol_[o];
		Neighs[cindex[o]]=Neighs_[o];
		Surface[cindex[o]]=Surface_[o];
	}
}

void VoronoiBinary::bPrintBody(ofstream & fout){
	bool have_clusters=!Clusters.empty();
	fout.write(as_byte(have_clusters),sizeof(have_clusters));
	fout.write(as_byte(nframe),sizeof(nframe));
	if(have_clusters){
		dmpVector(fout,atClusters);
		dmpVector(fout,VolClusters);
		dmpVector(fout,Clusters);
		dmpVector(fout,AreaClusters);
	}
	vector<double> Vol_(cindex.size());
	vector<vector<int>> Neighs_(cindex.size());
	vector<vector<double>> Surface_(cindex.size());
	for(size_t o{0};o<cindex.size();o++){
		Vol_[o]=Vol[cindex[o]];
		Neighs_[o]=Neighs[cindex[o]];
		Surface_[o]=Surface[cindex[o]];
	}
	dmpVector(fout,Vol_);
	dmpVector(fout,Neighs_);
	dmpVector(fout,Surface_);
	nframe++;
}
void VoronoiBinary::bReadHeader(ifstream & fin){
	fin.seekg (0, fin.end);
	int length = fin.tellg();
	fin.seekg (0, fin.beg);
	try{
	  if(!length) throw string("\n Something is wrong here!! File length is zero. \n");
	}catch(const string & s){
	  cout << s <<endl;
	  Finale::Finalize::Final();
	}
	fin.read(as_byte(nresid),sizeof(nresid));
	fin.read(as_byte(nr),sizeof(nr));
	fin.read(as_byte(nc),sizeof(nc));
	rdVector(fin,SelectedResidues);
	types=vector<int>(nr);
	atTypes=vector<int>(nr);

	Rdii=vector<double>(nr);
	fin.read(as_byte(types[0]),sizeof(types[0])*nr);
	fin.read(as_byte(atTypes[0]),sizeof(atTypes[0])*nr);
	fin.read(as_byte(Rdii[0]),sizeof(Rdii[0])*nr);
	rdVector(fin,Residue);
	rdVector(fin,TypesName);
	rdVector(fin,typesResidueMask);
	rdVector(fin,cindex);
	rdVector(fin,CIndex);
};

void VoronoiBinary::bPrintHeader(ofstream & fout){

	fout.write(as_byte(nresid),sizeof(nresid));
	fout.write(as_byte(nr),sizeof(nr));
	fout.write(as_byte(nc),sizeof(nc));
	dmpVector(fout,SelectedResidues);
	fout.write(as_byte(types[0]),sizeof(types[0])*nr);
	fout.write(as_byte(atTypes[0]),sizeof(atTypes[0])*nr);
	fout.write(as_byte(Rdii[0]),sizeof(Rdii[0])*nr);
	dmpVector(fout,Residue);
	dmpVector(fout,TypesName);
	dmpVector(fout,typesResidueMask);
	dmpVector(fout,cindex);
	dmpVector(fout,CIndex);
	cout << CIndex.size()<< " "<<CIndex.size()<< " "<<Residue.size()<< " "<<TypesName.size()<< " "<<Rdii.size()<< " "<<endl;
}
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
