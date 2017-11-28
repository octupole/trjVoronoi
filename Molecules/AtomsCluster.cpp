/*
 * AtomsCluster.cpp
 *
 *  Created on: Nov 28, 2017
 *      Author: marchi
 */

#include "AtomsCluster.h"
template <typename T>
int AtomsCluster<T>::calls=0;
template <class T>
vector<DDvect<T>> * PBCvect<T>::nxyz=nullptr;


template <typename T>
void AtomsCluster<T>::SetupPercolate(){
}
template <typename T>
void AtomsCluster<T>::SetupPercolate(Topol_NS::Topol & myTop){
	auto Reference=myTop.gReferenceResidues();
	auto Index=myTop.gCIndex();
	auto atmss=myTop.getAtomName();
	auto resn=myTop.getResinfo();
	vector<vector<int>> MySel;
	for(size_t o{0};o<Reference.size();o++){
		MySel.push_back(Index[Reference[o]]);
	}
	Perco=new Percolation<T>(MySel,rd,resn,atmss);
}


template <typename T>
void AtomsCluster<T>::Percolate() {
	try{
		if(!Perco) throw string("Should initialize percolation. Abort.");
	} catch(const string & s){
		cout << s <<endl;
		exit(1);
	}
	vector<Dvect> v{x};
	Matrix co(Mt.getCO());
	Matrix oc(Mt.getOC());
	Perco->doContacts(v,co,oc);
	Perco->gCluster();
	Perco->Accumulate();
	auto Male=Perco->getCluster();

	cout << Male.size()<<endl;
}
template<typename T>
void AtomsCluster<T>::Reconstruct(const string & y, TopolMic & MyTop){
	try{
		if(!MyTop.CheckResidue(y)) throw string(" Reference residue " + y + " does not exist in the system. Abort");
	} catch(const string & s){
		cout << s << endl;
		exit(1);
	}
	if(!calls++) {PBCvect<T> dummy;}
	vector<Dvect> nxyz=PBCvect<T>::getVec();
	typedef map<string,vector<vector<int> > >  ResMap;
	typedef ResMap::iterator ResIter;

	ResMap & ResDef=MyTop.getDef();
	vector<vector<int> > & Res=ResDef[y];
	Matrix co=Mt.getCO();
	Matrix oc=Mt.getOC();
	double v[DIM];
	v[XX]=co[XX][XX];v[YY]=co[YY][YY];v[ZZ]=co[ZZ][ZZ];
	double MinCO=*min_element(v,v+3)*0.5;  // Half the smallest axis is the cutoff distance for atom shifts
	Dvect Xref;

	// Apply boundary condition on every molecule as the default boundaries n .xtc files are atomic

	for(ResIter it=ResDef.begin(); it != ResDef.end();it++){
		vector<vector<int> > & Res0=it->second;
		for(unsigned int o=0;o<Res0.size();o++){
			Xref=xa[Res0[o][0]];
			for(unsigned int p=1;p<Res0[o].size();p++) {
				Dvect x1=Xref-xa[Res0[o][p]];
				Dvect xc1=co*x1;
				if(xc1.Norm() > MinCO) {
					Dvect tmp=*min_element(nxyz.begin(),nxyz.end(),CellComp<T>(x1,co));
					for(int q=0; q < DIM;q++) xa[Res0[o][p]][q]=xa[Res0[o][p]][q]-tmp[q];
				}
			}
		}
	}
	// > Take the reference type molecules and eliminate boundary conditions --> begin

	// >>  Compute center of mass for each molecule


	vector<Dvect> xaa;
	for(unsigned int o=0; o< Res.size();o++){
		Dvect cm=0;
		for(unsigned int p=0; p< Res[o].size();p++){
			int pp=Res[o][p];
			cm+=xa[pp];
		}
		cm/=static_cast<double>(Res[o].size());
		xaa.push_back(cm);
	}


	// >> Obtain the center of mass for the cluster of the reference molecules

	Xref=xaa[0];
	for(unsigned n=1;n < xaa.size();n++){
		Dvect xb=xaa[n];
		Dvect tmp=*min_element(nxyz.begin(),nxyz.end(),CellComp<T>(Xref,xb,co));
		xaa[n]=xaa[n]-tmp;
	}

	Dvect cma=T{0.0};
	for(unsigned n=0;n < xaa.size();n++){
		cma+=xaa[n];
	}
	cma/=static_cast<double> (xaa.size());
	// >>Shift all system by cma

	for(int n=0;n<nr;n++)
		for(int m=0;m<DIM;m++)
			xa[n][m]=xa[n][m]-cma[m];


	//> Apply molecular boundary conditions on the entire system
	for(ResIter it=ResDef.begin(); it != ResDef.end();it++){
		vector<vector<int> > & Res0=it->second;
		for(unsigned int o=0;o<Res0.size();o++){
			Dvect xcm=T{0.0};
			for(unsigned int p=0;p< Res0[o].size();p++) {
				xcm+=xa[Res0[o][p]];
			}

			xcm/=static_cast<double>(Res0[o].size());

			for(int m=0;m<DIM;m++) xcm[m]=rint(xcm[m]);
			for(unsigned int p=0;p<Res0[o].size();p++)
				for(int m=0;m<DIM;m++)
					xa[Res0[o][p]][m]=xa[Res0[o][p]][m]-xcm[m];

		}
	}

		//> Obtain new Cartesian coordinates from reduced xa's
	for(int i=0;i<nr;i++){
		for(int o=0;o<DIM;o++){
			xa[i][o]=xa[i][o]+0.5;
			x[i][o]=Mt.getCO()[o][XX]*xa[i][XX]+Mt.getCO()[o][YY]*xa[i][YY]+Mt.getCO()[o][ZZ]*xa[i][ZZ];
		}
	}
}

template<typename T>
void AtomsCluster<T>::__ReconstructOneCluster(vector<bool> & atSolv){
	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();
	Matrix co=Mt.getCO();
	Matrix oc=Mt.getOC();
	Dvect xcmC{0};
	size_t u{0};
	for(size_t p=0;p<mCluster[0].size();p++){
		int n=mCluster[0][p];
		for(size_t i=0;i<mAtoms[n].size();i++){
			int ia=mAtoms[n][i];
			xcmC[XX]+=xa[ia][XX];
			xcmC[YY]+=xa[ia][YY];
			xcmC[ZZ]+=xa[ia][ZZ];
			u++;
		}
	}
	xcmC/=static_cast<double>(u);
	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			 xa[ia][o]-=xcmC[o]-HALF;
			 if(atSolv[ia]) {
				 xa[ia][o]-=rint(xa[ia][o]-HALF);
			 }
		}
	}
	//> Obtain new Cartesian coordinates from reduced xa's
	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			x[ia][o]=Mt.getCO()[o][XX]*xa[ia][XX]+Mt.getCO()[o][YY]*xa[ia][YY]+Mt.getCO()[o][ZZ]*xa[ia][ZZ];
		}
	}
}

template<typename T>
void AtomsCluster<T>::Reconstruct(Contacts<T> * con0){
	if(!calls++) {PBCvect<T> dummy;}
	vector<Dvect> nxyz=PBCvect<T>::getVec();
	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();

	Matrix co=Mt.getCO();
	Matrix oc=Mt.getOC();
	double v[DIM];
	v[XX]=co[XX][XX];v[YY]=co[YY][YY];v[ZZ]=co[ZZ][ZZ];
	double MinCO=*min_element(v,v+3)*0.5;  // Half the smallest axis is the cutoff distance for atom shifts
	Dvect Xref;

	// Reconstruct each one of the molecules
	vector<Dvect> xcm(mAtoms.size(),T{0.0});
	auto PBC=[](Dvect xa,Dvect xb){
		Dvect tmp{0};
		for(auto p=0;p<DIM;p++)
			tmp[p]=-rint(xa[p]-xb[p]);
		return tmp;};
	for(size_t o=0;o<nr;o++){
		for(size_t q=0;q<DIM;q++) xa[o][q]=xa[o][q]-rint(xa[o][q]);
	}
	for(size_t o=0;o<mAtoms.size();o++){

		Xref=xa[mAtoms[o][0]];
		for(size_t p=0;p<mAtoms[o].size();p++){
			int n=mAtoms[o][p];

			Dvect x1{Xref-xa[n]};
			Dvect tmp=*min_element(nxyz.begin(),nxyz.end(),CellComp<T>(x1,co));
			//			Dvect tmp{PBC(Xref,xa[n])};
			for(int q=0; q < DIM;q++) xa[n][q]=xa[n][q]-tmp[q];
			Xref=xa[n];
		}
		// Compute molecule center of mass xcm
		for(size_t p=0;p<mAtoms[o].size();p++){
			Dvect xx=xa[mAtoms[o][p]];
			xcm[o]+=xx;
		}
		xcm[o]/=static_cast<double>(mAtoms[o].size());
	}
	vector<bool> atSolv(nr,true);
	for(auto o=0;o<mAtoms.size();o++){
		for(auto p=0;p<mAtoms[o].size();p++){
			atSolv[mAtoms[o][p]]=false;
		}
	}

	//
	vector<Dvect> xcmC(mCluster.size());


	for(size_t o=0;o<mCluster.size();o++){
//		if(mCluster[o].size() < 5) continue;
		vector<Dvect> xcm0(mCluster[o].size());
		for(size_t p=0;p<mCluster[o].size();p++){
			xcm0[p]=xcm[mCluster[o][p]];
		}
		(*con0)(xcm0,co,oc);
		con0->Neighbors();
		vector<vector<int>> & nnl=con0->NNL();

// Start with a reference molecule p=0;

		for(size_t p=0;p<mCluster[o].size();p++){
			size_t p0=mCluster[o][p];
			Xref=xcm0[p]; // Initial reference molecule is the first on the list
			for(size_t r=0;r<nnl[p].size();r++){
				size_t n=nnl[p][r];
				size_t n0=mCluster[o][n];
				Dvect x1{Xref-xcm0[n]};
				Dvect tmp=*min_element(nxyz.begin(),nxyz.end(),CellComp<T>(x1,co));
				//				Dvect tmp{PBC(Xref,xcm0[n])};
				if(p>n) continue;
				for(int q=0; q < DIM;q++) xcm0[n][q]=xcm0[n][q]-tmp[q];
				for(int q=0; q < DIM;q++) xcm[n0][q]=xcm[n0][q]-tmp[q];
				for(size_t i=0;i<mAtoms[n0].size();i++){
				  int ia=mAtoms[n0][i];
				  for(int q=0; q < DIM;q++) xa[ia][q]=xa[ia][q]-tmp[q];
				}

			}
		}


		// Compute the center of mass of the cluster
		xcmC[o]=0.0;
		for(size_t p=0;p<mCluster[o].size();p++){
			xcmC[o]+=xcm0[p];
		}
		xcmC[o]/=static_cast<double>(mCluster[o].size());
	}
	if(mCluster.size() == 1) return __ReconstructOneCluster(atSolv);

	//> Translate the cell to get all images inside the cell.
	Dvect Translation{__FindCell(mCluster,mAtoms)};

	//> Translate all coordinates as the center of the super cluster must be the center of the system

	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			xa[ia][o]+=Translation[o];
			if(atSolv[ia]) xa[ia][o]-=rint(xa[ia][o]-HALF);
		}
	}
	for(size_t o=0;o<mAtoms.size();o++){
		// Compute molecule center of mass xcm
		xcm[o]=0.0;
		for(size_t p=0;p<mAtoms[o].size();p++){
			Dvect xx=xa[mAtoms[o][p]];
			xcm[o]+=xx;
		}
		xcm[o]/=static_cast<double>(mAtoms[o].size());
	}
	for(size_t o=0;o<mCluster.size();o++){
		Dvect xcmC{T{0.0}};
		for(size_t p=0;p<mCluster[o].size();p++){
			xcmC+=xcm[mCluster[o][p]];
		}
		xcmC/=static_cast<double>(mCluster[o].size());
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t i=0;i<mAtoms[n].size();i++){
				int ia=mAtoms[n][i];
				xa[ia][XX]=xa[ia][XX]-rint(xcmC[XX]-HALF);
				xa[ia][YY]=xa[ia][YY]-rint(xcmC[YY]-HALF);
				xa[ia][ZZ]=xa[ia][ZZ]-rint(xcmC[ZZ]-HALF);
			}

		}
	}

	//> Obtain new Cartesian coordinates from reduced xa's
	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			x[ia][o]=Mt.getCO()[o][XX]*xa[ia][XX]+Mt.getCO()[o][YY]*xa[ia][YY]+Mt.getCO()[o][ZZ]*xa[ia][ZZ];
		}
	}
}
template <typename T>
DDvect<T> AtomsCluster<T>::__FindCell(const vector<vector<int> > & mCluster, const vector<vector<int> > & mAtoms){
	vector<Dvect> x(nr);
	Matrix CO(Mt.getCO());
	for(auto o=0;o<nr;o++){
		for(auto n=0;n<DIM;n++){
			x[o][n]=xa[o][n];
		}
	}
	vector<double> xcm(mAtoms.size(),0.0);
	Dvect myReturn{T{0.0}};
	double dx=0.20;
	for(auto n=0;n<DIM;n++){
		size_t Nn=CO[n][n]/dx;
		size_t M{0};
		double ddx=1.0/Nn;
		vector<int> myPBC(Nn);
		while(M < Nn){
			size_t bCount{0};
			for(auto o=0;o<nr;o++){
				x[o][n]+=ddx;
			}
			for(size_t o=0;o<mAtoms.size();o++){
				// Compute molecule center of mass xcm
				xcm[o]=0.0;
				for(size_t p=0;p<mAtoms[o].size();p++){
					double xx=x[mAtoms[o][p]][n];
					xcm[o]+=xx;
				}
				xcm[o]/=static_cast<double>(mAtoms[o].size());
				xcm[o]-=rint(xcm[o]-0.5);
			}

			for(size_t o=0;o<mCluster.size();o++){
				int nn=mCluster[o][0];
				double Xref=xcm[nn];
				for(size_t p=1;p<mCluster[o].size();p++){
					int nn=mCluster[o][p];
					double PBC=rint(Xref-xcm[nn]);
					if(PBC) bCount++;
				}
			}
			myPBC[M]=bCount;
			M++;
		}
		auto result = std::minmax_element(myPBC.begin(),myPBC.end());
		auto it1=std::find_if(myPBC.begin(),myPBC.end(),[result](int q){return q==*result.first;});
		auto it2=std::find_if(it1,myPBC.end(),[result](int q){return q!=*result.first;});
		auto j=std::distance(myPBC.begin(),it1)+std::distance(it1,it2)/2;
		myReturn[n]=j*ddx;
	}
	return 	myReturn;
}


template <typename T>
AtomsCluster<T>::~AtomsCluster() {
	delete Perco;
}
template class AtomsCluster<float>;
template class AtomsCluster<double>;
