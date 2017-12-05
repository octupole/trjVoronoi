/*
 * Topol.cpp
 *
 *  Created on: May 12, 2012
 *      Author: marchi
 */

#include "TopolMic.h"

namespace Topol_NS {
TopolMic::TopolMic():Topol::Topol(){

}
TopolMic::TopolMic(TopolPDB & data, bool x): Topol::Topol(data,x) {
}
void TopolMic::ExtractInfoMic(TopolPDB & data,bool bRD){
	//Topol::ExtractInfo(data,bRD);
	int nat=0;
	MyResidue mytmp;
	map<string,map<int,vector<int> > > Defs,DefsNH;

	for(size_t i=0;i<data.Size();i++){
		string sub1=data[i].atn; // Atom name
		string sub2=data[i].resn; // Residue name
		int nn=data[i].res;
		Defs[sub2][nn].push_back(nat);
		if(sub1.find_first_not_of("1234") == 1){
			if(sub1.compare(1,1,"H") != 0){
				DefsNH[sub2][nn].push_back(nat);
			}
		} else
			if(sub1.compare(0,1,"H") != 0){
				DefsNH[sub2][nn].push_back(nat);
			}

		nat++;
	}
	typedef map<string,map<int,vector<int> > > mappa;
	mappa::iterator itt=Defs.begin();
	for(;itt != Defs.end();itt++){
		map<int,vector<int> > & idx=itt->second;
		DefRes[itt->first]=vector<vector<int> >(idx.size());
		map<int,vector<int> >::iterator ipp=idx.begin();
		int ia=0;
		for(;ipp != idx.end() ;ipp++){

			DefRes[itt->first][ia]=ipp->second;
			ia++;
		}
	}
	itt=DefsNH.begin();
	for(;itt != DefsNH.end();itt++){

		map<int,vector<int> > & idx=itt->second;
		DefResNH[itt->first]=vector<vector<int> >(idx.size());
		map<int,vector<int> >::iterator ipp=idx.begin();
		int ia=0;
		for(;ipp != idx.end() ;ipp++){
			DefResNH[itt->first][ia]=ipp->second;
			ia++;
		}
	}
}


TopolMic::~TopolMic() {
	// TODO Auto-generated destructor stub
}
bool TopolMic::CheckResidue(const string & y){
	return DefRes.count(y) > 0;
}


} /* namespace Topol_NS */
