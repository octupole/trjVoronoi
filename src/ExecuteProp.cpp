/*
 * ExecuteVoronoi.cpp
 *
 *  Created on: Nov 16, 2017
 *      Author: marchi
 */

#include "ExecuteProp.h"
#include "Atoms.h"

namespace Properties {

template <typename T>
size_t ExecuteProp<T>::nnx=1;
template <typename T>
size_t ExecuteProp<T>::nny=1;
template <typename T>
size_t ExecuteProp<T>::nnz=1;
template <typename T>
Parallel::NewMPI * ExecuteProp<T>::CurrMPI{nullptr};
template <typename T>
void ExecuteProp<T>::__Print(ofstream & y){}

template <typename T>
ExecuteProp<T>::ExecuteProp(trj::TrjRead & MyIn) {
	__SetUp(MyIn);

	bOnce=MyIn.bbOnce();
	nnx=MyIn.gnnx();
	nny=MyIn.gnny();
	nnz=MyIn.gnnz();
	nstart=MyIn.gnstart();
	nend=MyIn.gnend();
	nskip=MyIn.gnskip();
	bTest=MyIn.bbTestVol();

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

}
template <typename T>
ExecuteProp<T>::ExecuteProp(trj::TrjRead & MyIn, Topol & Topology):
 	 ExecuteProp<T>::ExecuteProp(MyIn){
	Top=&Topology;
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

	Percolation<T>::setPercoCutoff(MyIn.gPercoCutoff());

	Clustering=true;
	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();
};
template <typename T>
void ExecuteProp<T>::operator()(Atoms<T> * atx){

	if(finx)
		__RunTrajectory(atx);
	else
		__RunPDB(atx);
}

template <typename T>
void ExecuteProp<T>::__RunTrajectory(Atoms<T> * atmx){

	myiterators::IteratorAtoms<T> iter_atm(atmx,finx,nstart,nend,nskip);

	while((++iter_atm).isReferenced()){
		stringstream ss;
		Atoms<T> * atmA=iter_atm();

		float ntime=atmA->getTime();
		int nClusters{0};
		atmA->setTopol(*Top);
		if(Clustering){

			static struct Once{
				Once(Atoms<T> * atmA, Topol_NS::Topol * myTop, bool whichJSON){
					atmA->SetupPercolate(*myTop,whichJSON);

				}
			} _Once(atmA, Top,JSONOutput);
			if(bOnce){
				static struct Once_p{int nClusters{0};Once_p(Atoms<T> *atmA){nClusters=atmA->Percolate();}} __Once_p(atmA);
				nClusters=__Once_p.nClusters;
			}else {
				nClusters=atmA->Percolate();
			}
		}
		switch(nClusters){
		case 0:
			break;
		case 1:
			ss << "    " <<fixed << setw(4) << nClusters<<" cluster  <-----";
			break;
		default:
			ss<< "    " << fixed << setw(4) << nClusters<<" clusters <-----";
		}

		atmA->Gyro();

		Comms->getStream() << atmA->getRg_i();
		Comms->getStream() << atmA->getComp();

		cout << fixed << setw(5) << "----> Time Step " << ntime << ss.str()<<"\n";
	}

	Comms->appendStreams();
	if(bDel) Comms->removeFiles();
	Comms->closeStream();   // close stream!!!
	CurrMPI->Barrier();
	CurrMPI->~NewMPI();
	cout << "\nProgram completed: Output data written to " + fileout << "\n\n";

}
template <typename T>
void ExecuteProp<T>::__RunPDB(Atoms<T> * atm){
	bool bTestVol{true};
	vector<string> data;
	fpdb->clear();
	fpdb->seekg(0);
	for(string str;getline(*fpdb,str);){
		data.push_back(str);
	}
	atm->pdb(data);
	float ntime=atm->getTime();
	atm->Gyro();

	Comms->getStream() << atm->getRg_i();
	Comms->getStream() << atm->getComp();
	Comms->appendStreams();
	if(bDel) Comms->removeFiles();
	CurrMPI->Barrier();
	CurrMPI->~NewMPI();
	cout << "\nProgram completed: Output data written to " + fileout << "\n\n";
}


template <typename T>
void ExecuteProp<T>::__SetUp(trj::TrjRead & MyIn){
	fileout=MyIn.gfileout();

	bDel=MyIn.bbDel();
	bHyd=MyIn.bbHyd();
	binOutput=MyIn.bbOutBin();
	JSONOutput=MyIn.bbOutJSON();
	fileout_bin=MyIn.gfileout_bin();

	fpdb=MyIn.gFpdb();
	finx=MyIn.gFinx();
	fin1x=MyIn.gFin1();
	foutx=MyIn.gFoutx();

	fidb=MyIn.gFidb();
	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();
	foutx=new ofstream();


	/*
	 * Define and dimension RhoSaxs and Saxs classes
	 */

	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();
	Clustering=true;
}

template <typename T>
ExecuteProp<T>::~ExecuteProp() {}
template <typename T>
ofstream & operator<<(ofstream & fout, ExecuteProp<T> & y)
{
		y.__Print(fout);
		return fout;
	}
template class ExecuteProp<float>;
template class ExecuteProp<double>;

} /* namespace Voronoi */
