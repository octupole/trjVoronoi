/*
 * Voronoi.cpp
 *
 *  Created on: Jun 20, 2011
 *      Author: marchi
 */

#include "Voronoi.h"

namespace Voro{
int Voronoi::nresid=0;
int Voronoi::nr=0;
int Voronoi::nc=0;
int Voronoi::nframe=0;
float Voronoi::time=0.0;
vector<string> Voronoi::label=vector<string>();





Voronoi::Voronoi() {}

Voronoi::Voronoi(Topol & myTop, bool bH){
	nr=myTop.Size();
	Vol=vector<double>(nr);
	Neighs=vector<vector<int>>(nr);
	Surface=vector<vector<double>>(nr);
	label=vector<string>(nr);
	types=vector<int>(nr);
	atTypes=vector<int>(nr);
	SelectedResidues=myTop.gReferenceResidues();


	label=myTop.getAtomName();
	vector<vector<int> > & Index=myTop.gCIndex();
	CIndex=vector<vector<int>>(Index.size());
	for(size_t o=0;o<Index.size();o++){
		for(size_t p=0;p<Index[o].size();p++){
			size_t i=Index[o][p];
			Vol[i]=0.0;
			if(!bH){
				string sub1=label[i];
				if(sub1.find_first_not_of("1234") == 1){
					if(sub1.compare(1,1,"H") != 0) {
						cindex.push_back(i);
						CIndex[o].push_back(i);
					}
				} else{
					if(sub1.compare(0,1,"H") != 0) {
						cindex.push_back(i);
						CIndex[o].push_back(i);
					}
				}
			} else {
				cindex.push_back(i);
				CIndex[o].push_back(i);
			}
		}
	}
	if(bH) {
		cout << std::left << std::fixed<< "\n              " <<std::setw(10) <<
 				" Voronoi's calculation includes hydrogens " << endl;
	} else{
		cout << std::left << std::fixed<< "\n              " <<std::setw(10) <<
 				" Voronoi's calculation does not includes hydrogens " << endl;
	}
 	cout <<"\n" << std::left << std::fixed<<"                  "<<std::setw(10) << " Atoms are  = " << "      "
 			<< std::setw(12) << cindex.size() << "\n" <<endl;

	atTypes=myTop.gatResType();
	nresid=myTop.ResSize();

	TypesName=myTop.getTypeNames();
	const vector<int> & tmp=myTop.getAtomTypeNo();
	for(int o=0;o<nr;o++)
		types[o]=tmp[o];
	RealResidue=myTop.getResind();
	Residue=myTop.getResinfo();
	Rdii=myTop.getrd();
	nc=myTop.getmaxt();
	area.Allocate(nresid,nc);
	Vols.Allocate(nresid);
	set<int> tmp1;
	for(string it: this->TypesName){
		tmp1.insert(ResidueTypes::find(it));
	}
	for(auto it:tmp1){
		this->typesResidueMask.push_back(it);
	}
}
void Voronoi::__extraInit(Topol & myTop, bool bH){
}
template <typename T>
void Voronoi::Start(float frame, Atoms<T> & atm){
		time=frame;
	const int nx=NNN,ny=NNN,nz=NNN;
	Metric<T> Mt=atm.getMt();
	MMatrix<T> CO=Mt.getCO();
	MMatrix<T> OC=Mt.getOC();
	T VolTmp=Mt.getVol();
	Matrix baba;
	for(int o{0};o<DIM;o++)
		for(int p{0};p<DIM;p++){
			co[o][p]=static_cast<double>(CO[o][p]);
			oc[o][p]=static_cast<double>(OC[o][p]);
		}
	VolCell=static_cast<double>(VolTmp);
	double bx=co[0][0], bxy=co[0][1], by=co[1][1],bxz=co[0][2],byz=co[1][2],bz=co[2][2];

	vector<double> vertx;
	if(Mycon) delete Mycon;
	if(porder) delete porder;
	Mycon=new container_periodic_poly(bx,bxy,by,bxz,byz,bz,nx,ny,nz,8);
	porder=new particle_order;
	for(unsigned int o=0;o<cindex.size();o++){
		double x=atm[cindex[o]][XX];
		double y=atm[cindex[o]][YY];
		double z=atm[cindex[o]][ZZ];
		double r=Rdii[cindex[o]];
		Mycon->put(*porder,o,x,y,z,r);
	}
	Percolation<T> * myPerco=atm.gPerco();
	if(myPerco){
		Clusters.clear();
		atClusters.clear();
		VolClusters.clear();
		AreaClusters.clear();
		SurfaceClusters.clear();
		listcon mAtoms=myPerco->getAtoms();
		listcon mClusters=myPerco->getCluster();
		Clusters=listcon(mClusters.size());
		atClusters=vector<int>(nr,-1);
		for(size_t o{0};o<mClusters.size();o++)
			for(size_t p{0};p<mClusters[o].size();p++){
				int i0{mClusters[o][p]};
				for(size_t n{0};n<CIndex[i0].size();n++){
					int i{CIndex[i0][n]};
					Clusters[o].push_back(i);
					atClusters[i]=o;
				}
			}
		VolClusters=vector<double>(Clusters.size(),0.0);
		AreaClusters=vector<double>(Clusters.size(),0.0);
		vector<double> vv(nc,0.0);
		SurfaceClusters=vector<vector<double>>(Clusters.size(),vv);
	}
}
void Voronoi::testVol(){
	double VorVol=0.0;
	for(unsigned int n=0;n<cindex.size();n++)
		VorVol+=Vol[cindex[n]];
	cout << setw(10) << setprecision(2) << scientific
			<< "Volume error is " << 1000.0*(VorVol-VolCell)
			<< " A^3 over " << setprecision(4) << fixed << 1000.0*VolCell << " A^3 "<< endl;

}


void Voronoi::gather(vector<int> & it){
	vector<int> tmp=it;
	for(unsigned int o=0;o<it.size();o++) {
		tmp[o]=cindex[it[o]];
	}
	it=tmp;
	tmp.clear();
}

void Voronoi::getData(){
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
	for(int o=0;o<nresid;o++) {
		vector<int> & cindex=CIndex[o];
		double sum_v=0.0;
		vector<tArea> cdx(cindex.begin(),cindex.end());
		sort(cdx.begin(),cdx.end(),tAcomp());
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
				area[o][getTypes(it->n)]+=it->a;
		}
		Vols[o]=sum_v;
	}

}
void Voronoi::WriteIt(std::ofstream & fout){
	fout << "######>> At step No. " << setw(10) << setprecision(2) << fixed<< time << endl;
	int po=0;
	for(int o=0;o<nresid;o++) {
		if(!VoronoiSetter::bPrintVols)continue;
		if(getTypesRes(o) != VoronoiSetter::pGroup && VoronoiSetter::pGroup != -1) continue;
		double a=Vols[o]*1000.0;
		string l=Residue[o];
		int rres=RealResidue[o]+1;
		(po%5)?fout << setw(10) << setprecision(4) << fixed<<  a << ' ' << setw(4) << l << ' ' << setw(5) << rres << ' ':
		  fout << endl << "%$VolRes " << setw(10) << setprecision(4) << fixed << a << ' ' << setw(4) << l  <<' ' << setw(5) << rres << ' ';
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
		int o_type=getTypesRes(o);
		VolSel[o_type]+=Vols[o]*1000.0;
		for(int p=0;p<nc;p++) {
			double a=area[o][p]*100.0;
			interface[o_type][p]+=a;
		}
	}

	if(VoronoiSetter::bPrintAreas ){
		fout << "# Area format:        ";
		for(int o=0;o<nc;o++) fout << setw(1) << fixed << " areatype(" << o << ") ";
		fout << endl;
		fout  << "%$AreaRes " ;
		int po=0;
		for(int o=0;o<nresid ;o++){
			string l=Residue[o];
			int o_type=getTypesRes(o);
			if(VoronoiSetter::pGroup != -1 && o_type != VoronoiSetter::pGroup) continue;
			fout  << right << setw(5)<< l << " " << setw(5) << fixed << RealResidue[o]+1 << ' ' << setw(3) << fixed << o_type << ' ' ;
			for(int p=0;p<nc;p++) {
				double a=area[o][p]*100.0;
				fout << setw(10) << setprecision(5) << fixed << a << ' ';
			}
			if(o+1==nresid) fout << endl;
			else if(!((po+1)%2)) fout << endl << "%$AreaRes " ;
			po++;
		}
	}

	fout << endl;
	fout << "#  Residue Types: ";
	for(int o=0;o<nc;o++) fout << setw(1) << fixed << " " << TypesName[o] << "  ["<<  o << "],  ";
	fout << endl;

	fout << "# Volume of selection " << endl <<"#      ";
	for(int o=0;o<nc;o++) fout << setw(1) << fixed << " Volumetype(" << o << ")  ";

	fout << endl;
	fout << "%$TotVol ";
	for(int o=0;o<nc;o++){
			fout << setw(15) << setprecision(4) << fixed << VolSel[o];
	}
	fout << endl;

	fout << "# Interface area of selection " << endl << "#             ";
	for(int o=0;o<nc;o++) fout << setw(1) << fixed << "  areatype(" << o << ")";
	fout << endl;
	for(int o=0;o<nc;o++){
		fout << "%$AreaTot " << setw(3) << o ;
		for(int p=0;p<o;p++)
			fout <<setw(13) << setprecision(3) << fixed << ' ';

		for(int p=o;p<nc;p++)
			fout << setw(13) << interface[o][p] ;
		fout << endl;
	}


}
template <typename T>
void Voronoi::dmpVector(ofstream & f, vector<T> & v){
	size_t ntmp=v.size();
	f.write(as_byte(ntmp),sizeof(ntmp));
	f.write(as_byte(v[0]),sizeof(v[0])*ntmp);
};
template <typename T>
void Voronoi::dmpVector(ofstream & f, vector<vector<T>> & v){
	size_t ntmp=v.size();
	f.write(as_byte(ntmp),sizeof(ntmp));
	size_t ntmp0;
	for(size_t o{0};o<ntmp;o++){
		ntmp0=v[o].size();
		f.write(as_byte(ntmp0),sizeof(ntmp0));
		f.write(as_byte(v[o][0]),sizeof(v[o][0])*ntmp0);
	}
};
void Voronoi::dmpVector(ofstream & f, vector<string> & v){
	string tmp_s;
	for(auto it: v)
		tmp_s+=it+" ";
	size_t ntmp=tmp_s.size();
	f.write(as_byte(ntmp),sizeof(ntmp));
	f.write(as_byte(tmp_s[0]),sizeof(tmp_s[0])*ntmp);

}
template <typename T>
void Voronoi::rdVector(ifstream & f,vector<T> & v){
	size_t ntmp;
	f.read(as_byte(ntmp),sizeof(ntmp));
	v.clear();
	v=vector<T>(ntmp);
	f.read(as_byte(v[0]),sizeof(v[0])*ntmp);
}
template <typename T>
void Voronoi::rdVector(ifstream & f,vector<vector<T>> & v){
	size_t ntmp;
	f.read(as_byte(ntmp),sizeof(ntmp));
	v.clear();
	v=vector<vector<T>>(ntmp);
	size_t ntmp0;
	for(size_t o{0};o<ntmp;o++){
		f.read(as_byte(ntmp0),sizeof(ntmp0));
		v[o]=vector<T>(ntmp0);
		f.read(as_byte(v[o][0]),sizeof(v[o][0])*ntmp0);
	}
}
void Voronoi::rdVector(ifstream & f,vector<string> & v){
	size_t ntmp;
	v.clear();
	f.read(as_byte(ntmp),sizeof(ntmp));
	string tmp_s(ntmp,' ');
	f.read(as_byte(tmp_s[0]),sizeof(tmp_s[0])*ntmp);
	v=split(tmp_s);
}
void Voronoi::bReadBody(ifstream & fin){
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

void Voronoi::bPrintBody(ofstream & fout){
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
void Voronoi::bReadHeader(ifstream & fin){
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
	types=vector<int>(nc);
	atTypes=vector<int>(nc);
	Rdii=vector<double>(nc);
	RealResidue=vector<int>(nresid);
	fin.read(as_byte(types[0]),sizeof(types)*nc);
	fin.read(as_byte(atTypes[0]),sizeof(atTypes)*nc);
	fin.read(as_byte(Rdii[0]),sizeof(Rdii)*nc);
	fin.read(as_byte(RealResidue[0]),sizeof(RealResidue[0])*nresid);
	rdVector(fin,Residue);
	rdVector(fin,TypesName);
	rdVector(fin,typesResidueMask);
	rdVector(fin,cindex);
	rdVector(fin,CIndex);
};

void Voronoi::bPrintHeader(ofstream & fout){
	fout.write(as_byte(nresid),sizeof(nresid));
	fout.write(as_byte(nr),sizeof(nr));
	fout.write(as_byte(nc),sizeof(nc));
	dmpVector(fout,SelectedResidues);
	fout.write(as_byte(types[0]),sizeof(types)*nc);
	fout.write(as_byte(atTypes[0]),sizeof(atTypes)*nc);
	fout.write(as_byte(Rdii[0]),sizeof(Rdii)*nc);
	fout.write(as_byte(RealResidue[0]),sizeof(RealResidue[0])*nresid);
	dmpVector(fout,Residue);
	dmpVector(fout,TypesName);
	dmpVector(fout,typesResidueMask);
	dmpVector(fout,cindex);
	dmpVector(fout,CIndex);
}

template void Voronoi::Start(float, Atoms<double> &);
template void Voronoi::Start(float, Atoms<float> &);
//template void Voronoi::dmpVector(ofstream &,vector<int>&);
//template void Voronoi::dmpVector(ofstream &,vector<double>&);
//template void Voronoi::rdVector(ifstream &,vector<int>&);
//template void Voronoi::rdVector(ifstream &,vector<double>&);
}
