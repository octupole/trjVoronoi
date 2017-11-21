/*
 * TopolPDB.h
 *
 *  Created on: Aug 24, 2013
 *      Author: marchi
 */

#ifndef TOPOLPDB_H_
#define TOPOLPDB_H_

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cctype>
#include "MyUtilClass.h"
#include "Split.h"
#include "ResidueTypes.h"

namespace Topol_NS {
using namespace DVECT;

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::map;
using std::stringstream;
using std::copy;
using Dvect=DDvect<double>;

struct PDBdata{
	int at, res;
	string atn,resn,domn,type;
	Dvect x;
	float domd;
	string first, last;
};

typedef map<string,map<string,vector<string> > > Mymapmap;
typedef map<string,vector<string > > Mymap;
class TopolPDB {
	struct op_ResNo{
		int n;
		int oldno;
		op_ResNo(): n(0), oldno(-1){};
		string operator()(const string & y){
			string out=y;
			int currno;
			std::stringstream(out.substr(22,4))>> currno; // Residue number
			if(oldno != currno) n++;
			oldno=currno;
			std::stringstream ss;ss<<std::fixed << std::setw(5) << n;
			out.replace(22,4,ss.str());
			return out;
		}
		int getSize() const {return oldno;};
	};
	vector<PDBdata> PDB;
	vector<PDBdata> SplitResidue(vector<string> &);
	int offset;
	int offsetRes;
	vector<string> DetPolSegment;
public:
	TopolPDB(): offset(0), offsetRes(0){};
	void operator()(const vector<string> &);
	void sDetPolsegment(string,string);
	PDBdata & operator[](size_t n){return PDB[n];}
	size_t Size() const {return PDB.size();}
	virtual ~TopolPDB();
};

} /* namespace Topol_NS */
#endif /* TOPOLPDB_H_ */
