/*
 * VoronoiMicelles.cpp
 *
 *  Created on: Nov 19, 2017
 *      Author: marchi
 */

#include <VoronoiMicelles.h>
namespace Voro{
VoronoiMicelles::VoronoiMicelles(Topol & myTop, bool bH): Voronoi::Voronoi(myTop,bH){
	__extraInit(myTop,bH);
}
void VoronoiMicelles::__extraInit(Topol & myTop, bool bH){
	nc=ResidueTypes::Size();
	area.Deallocate();
	area.Allocate(nresid,nc);
}
void VoronoiMicelles::getData(){
	c_loop_order_periodic vl(*Mycon,*porder);
	voronoicell_neighbor c;

	double vol0;
	vector<int> nei;
	vector<double> area0;
	vector<double> norm0;
	for(unsigned int n=0;n<cindex.size();n++){
		Vol[cindex[n]]=0.0;
		Neighs[cindex[n]].clear();
		Surface[cindex[n]].clear();
	}

	int ia=0;


	if(vl.start())
		do {
			if(Mycon->compute_cell(c,vl)) {
				vol0=c.volume();
				c.neighbors(nei);
				c.face_areas(area0);
				gather(nei);
				Vol[cindex[ia]]=vol0;
				Neighs[cindex[ia]]=nei;
				Surface[cindex[ia]]=area0;
			}
			ia++;
		} while(vl.inc());
	nei.clear();
	area0.clear();


	area=0.0;
	for(int o=0;o<nresid;o++) {
		vector<int> & cindex=CIndex[o];
		double sum_v=0.0;
		vector<tArea> cdx(cindex.begin(),cindex.end());
		std::sort(cdx.begin(),cdx.end(),tAcomp());
		for(unsigned int ia=0;ia<cdx.size();ia++){
			int i=cdx[ia].n;
			vector<tArea> nei0,nei;
			sum_v+=Vol[i];
			for(unsigned int p=0;p<Neighs[i].size();p++)
				nei0.push_back(tArea(Neighs[i][p],Surface[i][p]));

			sort(nei0.begin(),nei0.end(),tAcomp());
			set_difference(nei0.begin(),nei0.end(),cdx.begin(),cdx.end(),back_inserter(nei),tAcomp());

			vector<tArea>::iterator it=nei.begin();
			for(;it != nei.end(); ++it)
				area[o][types[it->n]]+=it->a;
		}
		Vols[o]=sum_v;
	}

}
void VoronoiMicelles::WriteIt(std::ofstream & fout){
	static bool firsttime{true};
	if(firsttime){
		fout << "#  Defined Residue Types: "<<endl;
		fout << "#                 ";
		for(int o=0;o<nc;o++) {
			if(nc-o == 1)
				fout  <<ResidueTypes::getType(o)<< "("<< setw(1)<< o <<") ";
			else
				fout  <<ResidueTypes::getType(o)<< "("<< setw(1)<< o <<"), ";
		}
		fout << "#                 ";
		fout << endl;
		firsttime=false;
	}
	fout << endl;
	fout << "######>> At step No. " << setw(10) << setprecision(2) << fixed<< time << endl;
	int po=0;
	for(size_t o0=0;o0<this->SelectedResidues.size() ;o0++) {
		int o=this->SelectedResidues[o0];
		if(!VoronoiPrint::bPrintVols)continue;
		if(getTypesRes(o) != VoronoiPrint::pGroup && VoronoiPrint::pGroup != -1) continue;
		double a=Vols[o]*1000.0;
		string l=Residue[o];
		(po%4)?fout << setw(10) << setprecision(4) << fixed<<  a << ' ' << setw(4) << l << ' ' << setw(5) << o+1<< ' ':
		  fout << endl << "%$VolRes " << setw(10) << setprecision(4) << fixed << a << ' ' << setw(4) << l  <<' ' << setw(5) << o+1 << ' ';
		po++;
	}
	fout << endl;

	array2<double> interface;
	interface.Allocate(nc,nc);
	interface=0.0;
	array1<double> VolSel;
	VolSel.Allocate(nc);
	VolSel=0.0;
	for(int o=0;o<nresid ;o++){
		int o_type=ResidueTypes::find(Residue[o]);
		VolSel[o_type]+=Vols[o]*1000.0;
		for(int p=0;p<nc;p++) {
			double a=area[o][p]*100.0;
			interface[o_type][p]+=a;
		}
	}

	if(VoronoiPrint::bPrintAreas ){
		fout << "# Area format:             ";
		for(size_t o0=0;o0<this->typesResidueMask.size() ;o0++) {
			int o=this->typesResidueMask[o0];
			fout << setw(1) << fixed << setw(8)<<ResidueTypes::getType(o)<< "(" <<setw(1)<< o << ") ";
		}

		fout << endl;

		for(size_t o0=0;o0<SelectedResidues.size() ;o0++){
			fout  << "%$AreaRes " ;
			int o{SelectedResidues[o0]};
			string l=Residue[o];
			int o_type=ResidueTypes::find(l);
			if(VoronoiPrint::pGroup != -1 && o_type != VoronoiPrint::pGroup) continue;
			fout  << right << setw(5)<< l << " " << setw(5) << fixed << ResidueTypes::getType(o_type) << ' ' << setw(3) << fixed << o+1 << ' ' ;
			for(size_t p0=0;p0<this->typesResidueMask.size();p0++) {
				int p=this->typesResidueMask[p0];
				double a=area[o][p]*100.0;
				fout << setw(11) << setprecision(5) << fixed << a << ' ';
			}
			fout << endl;
		}
	}

	fout << endl;

	fout << "# Volume of selection " << endl <<"#        ";
	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		fout << fixed<<setw(12)<<right<<ResidueTypes::getType(o)<<"(" << setw(1)<< o << ")";
	}

	fout << endl;
	fout << "%$TotVol ";
	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		fout << setw(15) << setprecision(4) << fixed << VolSel[o];
	}
	fout << endl;

	fout << "# Interface area of selection " << endl << "#           ";
	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		fout << fixed<<setw(10)<<right<<ResidueTypes::getType(o)<<"(" << setw(1)<< o << ")";
	}
	fout << endl;
	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		fout << "#   ";
		fout <<setw(5) <<fixed << left<<ResidueTypes::getType(o)<<"(" << setw(1)<< o << ")";
		for(size_t p0{0};p0<o0;p0++) {
			fout <<setw(13)  << " ";
		}
		for(size_t p0{o0};p0<this->typesResidueMask.size();p0++) {
			int p=this->typesResidueMask[p0];
			fout << right<<setw(13) << interface[o][p] ;
		}
		fout << endl;
	}

}

VoronoiMicelles::~VoronoiMicelles() {
	// TODO Auto-generated destructor stub
}

}
