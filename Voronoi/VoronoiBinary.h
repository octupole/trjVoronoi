/*
 * VoronoiBinary.h
 *
 *  Created on: Dec 2, 2017
 *      Author: marchi
 */

#ifndef VORONOI_VORONOIBINARY_H_
#define VORONOI_VORONOIBINARY_H_

#include "VoronoiMicelles.h"

namespace Voro {

class VoronoiBinary: public VoronoiMicelles {
	void WriteIt(std::ofstream &);
	void ReadIt(std::ifstream &);
	void bPrintBody(ofstream &);
	void bPrintHeader(ofstream &);
	void bReadHeader(ifstream &);
	void bReadBody(ifstream &);
public:
	VoronoiBinary(ifstream &);
	VoronoiBinary(ofstream &,Topol &,bool);

	virtual ~VoronoiBinary();
};

} /* namespace Voro */

#endif /* VORONOI_VORONOIBINARY_H_ */
