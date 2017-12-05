/*
 * Topol.h
 *
 *  Created on: May 12, 2012
 *      Author: marchi
 */

#ifndef TOPOLMIC_H_
#define TOPOLMIC_H_
#include "Topol.h"
namespace Topol_NS {

class TopolMic: public Topol {
	vector<string> CurrentPDB;
	void ExtractInfoMic(TopolPDB & ,bool);
public:
	TopolMic();
	TopolMic(TopolPDB &, bool);
	vector<string> & getStringTypeNames(){return TypeNames;};
	virtual ~TopolMic();
	bool CheckResidue(const string & );

	vector<vector<int> > & ResIndx(const string & x){return DefRes[x];}
	vector<vector<int> > & ResIndxNH(const string & x){return DefResNH[x];}
	const vector<string> & getPDB() {return CurrentPDB;}
};

} /* namespace Topol_NS */
#endif /* TOPOLMIC_H_ */
