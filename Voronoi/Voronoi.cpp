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
void Voronoi::doVoro__(){
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

}

void Voronoi::getData(){
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

template <typename T>
void Voronoi::dmpVector(ofstream & f, vector<T> & v){
	size_t ntmp=v.size();
	f.write(as_byte(ntmp),sizeof(ntmp));
	if(ntmp) f.write(as_byte(v[0]),sizeof(v[0])*ntmp);
};
template <typename T>
void Voronoi::dmpVector(ofstream & f, vector<vector<T>> & v){
	size_t ntmp=v.size();
	f.write(as_byte(ntmp),sizeof(ntmp));
	if(!ntmp) return;
	size_t ntmp0;
	for(size_t o{0};o<ntmp;o++){
		ntmp0=v[o].size();
		f.write(as_byte(ntmp0),sizeof(ntmp0));
		if(ntmp0) f.write(as_byte(v[o][0]),sizeof(v[o][0])*ntmp0);
	}
};
void Voronoi::dmpVector(ofstream & f, vector<string> & v){
	string tmp_s;
	for(auto it: v)
		tmp_s+=it+" ";
	size_t ntmp=tmp_s.size();

	f.write(as_byte(ntmp),sizeof(ntmp));
	if(ntmp) f.write(as_byte(tmp_s[0]),sizeof(tmp_s[0])*ntmp);
}
template <typename T>
void Voronoi::rdVector(ifstream & f,vector<T> & v){
	size_t ntmp;
	f.read(as_byte(ntmp),sizeof(ntmp));
	if(!ntmp) return;
	v.clear();
	v=vector<T>(ntmp);
	f.read(as_byte(v[0]),sizeof(T)*ntmp);
}
template <typename T>
void Voronoi::rdVector(ifstream & f,vector<vector<T>> & v){
	size_t ntmp;
	f.read(as_byte(ntmp),sizeof(ntmp));
	v.clear();
	if(!ntmp) return;
	v=vector<vector<T>>(ntmp);
	size_t ntmp0;
	for(size_t o{0};o<ntmp;o++){
		f.read(as_byte(ntmp0),sizeof(ntmp0));
		if(ntmp0) {
			v[o]=vector<T>(ntmp0);
			f.read(as_byte(v[o][0]),sizeof(T)*ntmp0);
		}
	}
}
void Voronoi::rdVector(ifstream & f,vector<string> & v){
	size_t ntmp;
	v.clear();
	f.read(as_byte(ntmp),sizeof(ntmp));
	if(!ntmp) return;
	string tmp_s(ntmp,' ');
	f.read(as_byte(tmp_s[0]),sizeof(tmp_s[0])*(ntmp));
	v=split(tmp_s);
}

template void Voronoi::Start(float, Atoms<double> &);
template void Voronoi::Start(float, Atoms<float> &);
template void Voronoi::dmpVector(ofstream &,vector<int>&);
template void Voronoi::dmpVector(ofstream &,vector<vector<int>>&);
template void Voronoi::dmpVector(ofstream &,vector<double>&);
template void Voronoi::dmpVector(ofstream &,vector<vector<double>>&);
template void Voronoi::rdVector(ifstream &,vector<int>&);
template void Voronoi::rdVector(ifstream &,vector<vector<int>>&);
template void Voronoi::rdVector(ifstream &,vector<double>&);
template void Voronoi::rdVector(ifstream &,vector<vector<double>>&);
}
