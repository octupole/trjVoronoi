/*
 * ExecuteVoronoi.cpp
 *
 *  Created on: Nov 16, 2017
 *      Author: marchi
 */

#include "ExecuteVoronoi.h"
#include "Atoms.h"
#include "Voronoi.h"
#include "VoronoiMicelles.h"

namespace Voro {

template <typename T>
size_t ExecuteVoronoi<T>::nnx=1;
template <typename T>
size_t ExecuteVoronoi<T>::nny=1;
template <typename T>
size_t ExecuteVoronoi<T>::nnz=1;
template <typename T>
Parallel::NewMPI * ExecuteVoronoi<T>::CurrMPI{nullptr};
template <typename T>
void ExecuteVoronoi<T>::__Print(ostream & y){}

template <typename T>
ExecuteVoronoi<T>::ExecuteVoronoi(trj::TrjRead & MyIn, Topol & Topology): Top(&Topology){
	__SetUp(MyIn);
	bOnce=MyIn.bbOnce();
	nnx=MyIn.gnnx();
	nny=MyIn.gnny();
	nnz=MyIn.gnnz();
	nstart=MyIn.gnstart();
	nend=MyIn.gnend();
	nskip=MyIn.gnskip();
	bTest=MyIn.bbTestVol();
	try{
		if(nnx ==1 || nny == 1|| nnz == 1) throw string("Grid dimensions are not set!");
	}catch(const string & s){
		cout << s <<endl;
		Finale::Finalize::Final();
	}
	if(finx){
		size_t TotFrame=finx->gFrameStep();
		try{
			if(int(TotFrame) <nend) {
				stringstream ss0,ss1;
				ss0<< nend; ss1<< TotFrame;
				throw string("\n\nThe requested end-point \"")+ss0.str()
									+string("\" goes beyond the last frame of the \ntrajectory \""+ss1.str()+string("\"!!\n\n"));
			}
		}catch(const string & s){
			cout << s << endl;
			Finale::Finalize::Final();
		}
	}
	try{
		if(finx && nend-nstart+1 < nskip) throw string("\nNumber of selected steps is smaller than -skip parameter. Change and rerun.\n");
	}catch(const string & s) {cout << s <<endl;Finale::Finalize::Final();}

	if(CurrMPI->AmI_Parallel()){

		auto MyTaskNo=CurrMPI->Get_Rank();
		try{
			if(finx && (nend-nstart+1)/nskip < int(CurrMPI->Get_Size())) throw string(" The number of CPUs is larger than the number "
					"of the trajectory steps. Change and rerun. ");
		}catch(const string & s) {cout << s <<endl;Finale::Finalize::Final();}
		auto ntot=(nend-nstart+1)/nskip;
		auto ntask=CurrMPI->Get_Size();
		auto nChunk=(ntot/ntask)*nskip;
		nend=nstart+nChunk*(MyTaskNo+1)-1;
		nstart=nstart+nChunk*MyTaskNo;

	}
	vor=new VoronoiMicelles(Topology,bHyd);
	Percolation<T>::setPercoCutoff(MyIn.gPercoCutoff());
	Clustering=MyIn.bbClust();
};
template <typename T>
void ExecuteVoronoi<T>::operator()(Atoms<T> * atx){
	if(finx)
		__RunTrajectory(atx);
	else
		__RunPDB(atx);
}

template <typename T>
void ExecuteVoronoi<T>::__RunTrajectory(Atoms<T> * atmx){

	myiterators::IteratorAtoms<T> iter_atm(atmx,finx,nstart,nend,nskip);
	Contacts<T> * Con0;
	if(Rcut_in < 0) Rcut_in=15.0;
	Con0=new Contacts<T>(Rcut,Rcut_in);
	while((++iter_atm).isReferenced()){
		Atoms<T> * atmA=iter_atm();
		Con0->setR(Rcut_in,Rcut_in);

		float ntime=atmA->getTime();
		if(Clustering){
			atmA->setrd(*Top);
			static struct Once{Once(Atoms<T> * atmA){atmA->SetupPercolate();}} _Once(atmA);
			if(bOnce){
				static struct Once_p{Once_p(Atoms<T> *atmA){atmA->Percolate();}} __Once_p(atmA);
			}else atmA->Percolate();
			atmA->Reconstruct(Con0);
		}

		vor->Start(ntime,*atmA);
		vor->getData();

		Comms->getStream() << *vor;
		if(bTest) vor->testVol();
		cout << fixed << setw(5) << "----> Time Step " << ntime <<"\n";

	}

	Comms->appendStreams();
	if(bDel) Comms->removeFiles();
	CurrMPI->Barrier();
	CurrMPI->~NewMPI();
	cout << "\nProgram completed: Output data written to " + fileout << "\n\n";

}
template <typename T>
void ExecuteVoronoi<T>::__RunPDB(Atoms<T> * atm){
	bool bTestVol{true};
	vector<string> data;
	fpdb->clear();
	fpdb->seekg(0);
	for(string str;getline(*fpdb,str);){
		data.push_back(str);
	}
	atm->pdb(data);
	float ntime=atm->getTime();
	vor-> Start(ntime,*atm);
	vor->getData();
	Comms->getStream() << *vor;

	if(bTestVol) vor->testVol();
	Comms->appendStreams();
	if(bDel) Comms->removeFiles();
	CurrMPI->Barrier();
	CurrMPI->~NewMPI();
	cout << "\nProgram completed: Output data written to " + fileout << "\n\n";
}


template <typename T>
void ExecuteVoronoi<T>::__SetUp(trj::TrjRead & MyIn){
	fileout=MyIn.gfileout();

	bDel=MyIn.bbDel();
	bHyd=MyIn.bbHyd();
	fileout_bin=MyIn.gfileout_bin();

	fpdb=MyIn.gFpdb();
	finx=MyIn.gFinx();
	foutx=MyIn.gFoutx();

	fidb=MyIn.gFidb();
	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();
	foutx=new ofstream();


	ios::streampos len;
	HeaderTrj header;
// Read header of dcd file
	try{
		if(finx) {
			finx->seekg(0,"end");
			len=finx->tellg();
			finx->seekg(0,"beg");
			*finx>>header;
			try{
				if(!header.check(Top->Size())) {
					stringstream ss0,ss1;
					ss0<< Top->Size();
					ss1<<header.getNatoms();
					throw string("Number of atoms in the pdb ("+ss0.str()+") and trajectory files ("
							+ss1.str()+") does not match!");}
			}
			catch(const string & s){cout << s<<endl;Finale::Finalize::Final();}
			ofstream & fout=*foutx;

			CurrMPI->Barrier();

			Comms=new Parallel::FComms(CurrMPI,fout,fileout,nstart,nend,header.getNFR(),nskip);
			nstart=Comms->getStart();
			nend=Comms->getEnd();
		}else{
			ofstream & fout=*foutx;
			CurrMPI->Barrier();
			nstart=1;nend=1;nskip=1;
			if(CurrMPI->Get_Size() > 1) throw string("Cannot compute Voronoi from a PDB file in Parallel!");
			Comms=new Parallel::FComms(CurrMPI,fout,fileout,nstart,nend,header.getNFR(),nskip);
			nstart=Comms->getStart();
			nend=Comms->getEnd();
		}
	} catch(const string & s){cout << s<<endl;
	Finale::Finalize::Final();exit(1);}
	/*
	 * Define and dimension RhoSaxs and Saxs classes
	 */

	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();
	Clustering=MyIn.bbClust();
}

template <typename T>
ExecuteVoronoi<T>::~ExecuteVoronoi() {}

template class ExecuteVoronoi<float>;
template class ExecuteVoronoi<double>;

} /* namespace Voronoi */
