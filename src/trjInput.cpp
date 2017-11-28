/*
 * trjInput.cpp
 *
 *  Created on: May 22, 2012
 *      Author: marchi
 */

#include "trjInput.h"

namespace trj {
trjInput::trjInput(int ntot,char ** v) {
	vector<string> in;

	inmap["-o"]=in;
	inmap["-dcd"]=in;
	inmap["-xtc"]=in;
	inmap["-b"]=in;
	inmap["-e"]=in;
	inmap["-skip"]=in;
	inmap["-pdb"]=in;
	inmap["-pvol"]=in;
	inmap["-nohyd"]=in;
	inmap["-hyd"]=in;
	inmap["-nodel"]=in;
	inmap["-del"]=in;
	inmap["-rd"]=in;
	inmap["-nord"]=in;
	inmap["-test"]=in;
	inmap["-help"]=in;
	inmap["-parea"]=in;
	inmap["-pgroup"]=in;
	inmap["-surface"]=in;
	inmap["-protein"]=in;
	inmap["-water"]=in;
	inmap["-channels"]=in;
	inmap["-solute"]=in;
	inmap["-select"]=in;
	inmap["-detP"]=in;
	inmap["-detO"]=in;
	inmap["-det"]=in;
	inmap["-clust"]=in;

	map<string,vector<string> >::iterator it=inmap.begin();
	for(int n=0;it!=inmap.end();++it,n++){
		Usage[n]=" ";
	}
	Usage[0]="\t -dcd file.dcd \n";
	Usage[1]="\t -xtc file.xtc\n";
	Usage[2]="\t -pdb file.pdb \n";
	Usage[3]="\t -o fileout \n";
	Usage[4]="\t- b No. First Frame \n";
	Usage[5]="\t -e No. Last Frame\n";
	Usage[6]="\t -skip No. skipped frames \n";
	Usage[9]="\t -nohyd // Do not include hydrogens in Voronoi [default] \n";
	Usage[8]="\t -hyd // Include hydrogens in Voronoi\n";
	Usage[11]="\t -nodel // Do not delete temporary files in MPI runs \n";
	Usage[10]="\t -del // Delete temporary files in MPI runs [default]\n";
	Usage[12]="\t -rd // Use atomic radii in the Voronoi calculation [default]\n";
	Usage[13]="\t -nord // Do not use atomic radii in the Voronoi calculation \n";
	Usage[14]="\t -solute <string selection>\n"
			"\t\tDefine the residue of the solute to be used in the Voronoi analysis.\n";

	Usage[15]="\t -test // compute total volume from Voronoi or simulation cell parameters\n";
	Usage[16]="\t -help // write some on line help \n";
	Usage[17]="\t -pvol // Print volume \n";
	Usage[18]="\t -parea // Print area \n";
	Usage[19]="\t -det filename// Define the polar segment of a detergent residue \n"
			"\t\tThe file filename contains two strings:\n"
			"\t\t   i) the first string is the name of the detergent residue \n"
			"\t\t   ii)the second string contains the polar atoms of the residue\n"
			"\t\t*Important*: The polar atoms in the pdb and topology *MUST* all occur\n"
			"\t\tbefore the hydrophobic atoms.\n";
	Usage[20]="\t -detP <string>// Define a name for the polar segment of a detergent residue \n";
	Usage[21]="\t -detO <string>// Define a name for the hydrophobic segment of a detergent residue \n";
	Usage[22]="\t -shell [2]<int n>// Compute number and average volume of the n-shell water\n"
			"\t\tneighbours around solute \n";
	Usage[23]="\t -clust <float cutoff [0.0] in Angstroems>\n"
			"\t\tDo clustering by percolation of the solute. This is used to make sure that the solute is at the\n"
			"\t\tcenter of the cell. Doing clustering is expensive, it is always better to have the solute in the\n"
			"\t\tmiddle already at simulation time. If cutoff is zero, the default, the cutoff is chosen according\n"
			"\t\tmiddle to the Lennard-Jones sigma parameter\n";

	int n=1;
	string key;
	for(;n<ntot;){
		string tmp0(v[n]);
		if(v[n][0] =='-'){
			key.assign(v[n]);
			if(inmap.find(key) != inmap.end()){
				if(inmap[key].empty()) inmap[key].push_back(tmp0);
			} else{
				unknownC.push_back(key);
				inmap[key].push_back(tmp0);
			}
		}
		else{
			inmap[key].push_back(tmp0);
		}
		n++;
	}

}

vector<string> trjInput::getUsage(){
	vector<string> use;
	map<int,string>::iterator it=Usage.begin();
	for(;it!=Usage.end();++it){
		use.push_back(it->second);
	}
	return use;
}
trjInput::~trjInput() {
	// TODO Auto-generated destructor stub
}

} /* namespace trj */
