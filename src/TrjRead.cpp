/*
 * TrjRead.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: marchi
 */

#include "TrjRead.h"
namespace trj {
Parallel::NewMPI * TrjRead::CurrMPI=nullptr;

TrjRead::TrjRead(int nv,char ** v): trjInput::trjInput(nv,v) {
	// TODO Auto-generated constructor stub
	string comm0(v[0]);
	size_t mypos=comm0.find_last_of("/")+1;
	size_t length=comm0.size()-mypos;
	string command=comm0.substr(mypos,length);
	string errmsg;
	string Usage="Usage:\t"+ command + "\n";
	vector<string> use=this->getUsage();
	vector<string> SelRes;
	string Reference;
	for(unsigned int n=0;n<use.size();n++)
		Usage+=use[n];
	Usage+="\n\t Default values in square brackets []\n";
	try{
		if(nv == 1) throw Usage;
		if(int m=this->bTest().size()) {
			errmsg=" Command(s) not found: ";
			for(int n=0;n<m;n++)
				errmsg+=this->bTest()[n]+"  ";
			errmsg+="\n"+Usage;
			throw errmsg;
			}
		}
	catch(const string & s){
		cout << s << endl;
		Finale::Finalize::Final();
	}

}
void TrjRead::Input(){
	ifstream ftest,fdefdomain;
#ifdef __INTEL  //  set the conversion factors by hand as until intel c++ compiler version 16 user-defined literals are not allowed
	struct dummy{
		double val;
	};
	dummy e0_conversion{0.694461547828659};
	ConvFactor=55.7036562713774;
#else
	auto bimbo=1.0_el/(1.0_nm*1.0_nm);

	auto fluct_en=(1.0/(4.0*M_PI))*bimbo*1.0_el*1.0_nm/eps0;
	auto Factor=fluct_en/(kt300/avogad) ;
	ConvFactor=Factor.val;

	auto FieldConv=1.0_Volt/1.0_nm;
	auto e0_conversion=4.0*M_PI*FieldConv*eps0/bimbo;
#endif
	int pGroup{-1}; // From 0 to ngrps, ngrps no. of selected groups
	bool bPrintVols{true};
	bool bPrintAreas{false};

	try{
		if(!inmap["-dcd"].empty()) {
			if(inmap["-dcd"].size() < 2) throw string("\n Filename expected for " + inmap["-dcd"][0] + " option \n");
			if(inmap["-dcd"].size() > 2) throw string("\n More than one entry for " + inmap["-dcd"][0] + " option \n");
			filein=inmap["-dcd"][1];
			ftest.open(filein.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open " + filein + "!!\n");
			ftest.close();

			finx=new FstreamF(filein);
			inputfile=true;
		}
		if(!inmap["-xtc"].empty()) {
			if(inmap["-xtc"].size() < 2) throw string("\n filename expected for " + inmap["-xtc"][0] + " option\n ");
			if(inmap["-xtc"].size() > 2) throw string("\n More than one entry for " + inmap["-xtc"][0] + " option \n");
			filein=inmap["-xtc"][1];
			ftest.open(filein.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open " + filein + "!!\n");
			ftest.close();
			finx=new FstreamC(filein,"rb");
			inputfile=true;
		}
		if(!inmap["-pdb"].empty()) {
			if(inmap["-pdb"].size() < 2) throw string("\n filename expected for " + inmap["-pdb"][0] + " option \n");
			if(inmap["-pdb"].size() > 2) throw string("\n More than one entry for " + inmap["-pdb"][0] + " option \n");
			filepdb=inmap["-pdb"][1];
			fpdb.open(filepdb.c_str(),ios::in);
			if(!fpdb) throw string("\n Cannot open " + filepdb + "!!\n");
		}
		if(!inmap["-select"].empty()) {
			if(inmap["-select"].size() != 2) throw string(" String of selected residues needed for " + inmap["-select"][0] + " option ");
			string selection=inmap["-select"][1];
			stringstream iss(selection);
			copy(istream_iterator<string>(iss),
					istream_iterator<string>(),
					back_inserter<vector<string> >(SelRes));
		}
		if(!inmap["-solute"].empty()) {
			if(inmap["-solute"].size() != 2) throw string(" Reference residue for Micelle centering needed " + inmap["-solute"][0] + " option ");
			string Ref=inmap["-solute"][1];
			stringstream iss(Ref);
			copy(istream_iterator<string>(iss),
					istream_iterator<string>(),
					back_inserter<vector<string> >(Reference));
		}
		if(!inmap["-det"].empty()) {
			if(inmap["-det"].size() != 2) throw string("A filename is needed for " + inmap["-det"][0] + " option ");
			string filedet=inmap["-det"][1];
			ifstream fdet(filedet.c_str(),ios::in);
			if(!fdet) throw string("\n Cannot open " + filedet + "!!\n");
			vector<string> data;
			for(string str;getline(fdet,str);){
				data.push_back(str);
			}
			DetResidue=data[0];
			PolarAtoms=data[1];
		}

		if(!inmap["-o"].empty()) {
			if(inmap["-o"].size() < 2) throw string("\n filename expected for " + inmap["-o"][0] + " option \n");
			if(inmap["-o"].size() > 2) throw string("\n More than one entry for " + inmap["-o"][0] + " option \n");
			fileout=inmap["-o"][1];
		}
		if(!inmap["-b"].empty()) {
			if(inmap["-b"].size() != 2) throw string(" Number of first frame needed for " + inmap["-b"][0] + " option ");
			stringstream(inmap["-b"][1])>> nstart;
		}
		if(!inmap["-e"].empty()) {
			if(inmap["-e"].size() != 2) throw string(" Number of last frame needed for " + inmap["-e"][0] + " option ");
			stringstream(inmap["-e"][1])>> nend;
		}
		if(!inmap["-skip"].empty()) {
			if(inmap["-skip"].size() != 2) throw string(" Number of skipped frames needed for " + inmap["-skip"][0] + " option ");
			stringstream(inmap["-skip"][1])>> nskip;
		}
		if(!inmap["-nohyd"].empty()) {
			if(inmap["-nohyd"].size() != 1) throw string(" No argument for " + inmap["-nohyd"][0] + " option ");
				bHyd=false;
		}
		if(!inmap["-hyd"].empty()) {
			if(inmap["-hyd"].size() != 1) throw string(" No argument for " + inmap["-hyd"][0] + " option ");
				bHyd=true;
		}
		if(!inmap["-del"].empty()) {
			if(inmap["-del"].size() != 1) throw string(" No argument for " + inmap["-del"][0] + " option ");
				bDel=false;
		}
		if(!inmap["-rcut"].empty()) {
			if(inmap["-rcut"].size() != 2) throw string(" Linked list cutoff needed for  " + inmap["-rcut"][0] + " option ");
			stringstream(inmap["-rcut"][1])>> Rcut;
			Rcut*=unit_nm;
		}
		if(!inmap["-nord"].empty()) {
			if(inmap["-nord"].size() != 1) throw string(" No argument for " + inmap["-nord"][0] + " option ");
				bIsrd=false;
		}
		if(!inmap["-pvol"].empty()) {
			if(inmap["-pvol"].size() != 1) throw string(" No argument for " + inmap["-pvol"][0] + " option ");
				bPrintVols=true;
		}
		if(!inmap["-pgroup"].empty()) {
			if(inmap["-pgroup"].size() != 2) throw string(" Group number needed for " + inmap["-pgroup"][0] + " option ");
			stringstream(inmap["-pgroup"][1])>> pGroup;
		}
		if(!inmap["-parea"].empty()) {
			if(inmap["-parea"].size() != 1) throw string(" No argument for " + inmap["-parea"][0] + " option ");
				bPrintAreas=true;
		}
		if(!inmap["-detP"].empty()) {
			if(inmap["-detP"].size() < 2) throw string("\n String expected for " + inmap["-detP"][0] + " option\n ");
			if(inmap["-detP"].size() > 2) throw string("\n More than one entry for " + inmap["-detP"][0] + " option \n");
			string myDetP=inmap["-detP"][1];
			Topol_NS::ResidueTypes::addDetPolar(myDetP);
		}
		if(!inmap["-detO"].empty()) {
			if(inmap["-detO"].size() < 2) throw string("\n String expected for " + inmap["-detO"][0] + " option\n ");
			if(inmap["-detO"].size() > 2) throw string("\n More than one entry for " + inmap["-detO"][0] + " option \n");
			string myDetO=inmap["-detO"][1];
			Topol_NS::ResidueTypes::addDetOil(myDetO);
		}
		if(!inmap["-test"].empty()) {
			if(inmap["-test"].size() != 1) throw string(" No argument for " + inmap["-test"][0] + " option ");
				bTestVol=true;
		}

		if(!inmap["-grid"].empty()) {
			int Nm=inmap["-grid"].size();
			try{
				if(Nm == 2){
					stringstream(inmap["-grid"][1])>> nnx;
					nny=nnx; nnz=nnx;
				} else if(Nm == 4){
					stringstream(inmap["-grid"][1])>> nnx;
					stringstream(inmap["-grid"][2])>> nny;
					stringstream(inmap["-grid"][3])>> nnz;
				} else throw string("\n Warning: No Grid points provided " + inmap["-grid"][0] + " option. "
						+"Hope this is correct \n");
			}catch(const string & s){
				cout << s <<endl;
			}
		}
		if(!inmap["-clust"].empty()) {
			if(inmap["-clust"].size() > 2) throw string("\n Only one optional argument is allowed for " + inmap["-clust"][0] + "\n");
			bClust=true;
			if(inmap["-clust"].size() == 2) {
				stringstream(inmap["-clust"][1])>> PercoCutoff;
				PercoCutoff/=10.0;
			}
		}
		if(!inmap["-once"].empty()) {
			if(inmap["-once"].size() != 1) throw string("\n No argument for " + inmap["-once"][0] + " option \n");
			bOnce=true;
		}
	}catch(const string & s){
		cout << s << endl;
		Finale::Finalize::Final();
	}

	string tmp1=fileout.substr(0,fileout.find_first_of("."));
	string tmp2=".dat";


	gFinx=finx;
	gFin2x=fin2x;
	gFout_xtcx=fout_xtcx;
	gFout_pdbx=fout_pdbx;
	gFoutx=foutx;

	gFpdb=&fpdb;
	VoronoiPrint::bPrintAreas=bPrintAreas;
    VoronoiPrint::bPrintVols=bPrintVols;
    VoronoiPrint::pGroup=pGroup;

	try{
		if(!Reference.empty() && !SelRes.empty()) {
			vector<string> ss;
			std::sort(SelRes.begin(),SelRes.end());std::sort(Reference.begin(),Reference.end());
			std::set_intersection(SelRes.begin(),SelRes.end(),Reference.begin(),Reference.end(),back_inserter(ss));
			if(ss.empty()) throw string("Residue in -solute and -select must intersect");
		}
		if(finx && nstart > nend) throw string("Initial trajectory step is larger than the set final step. Change and rerun. ");
	} catch(const string & s){
		cout << s <<endl;
		Finale::Finalize::Final();
	}
	Mycut/=unit_nm;
	Myd/=unit_nm;
}
TrjRead::~TrjRead() {
	// TODO Auto-generated destructor stub
}

} /* namespace trj */
