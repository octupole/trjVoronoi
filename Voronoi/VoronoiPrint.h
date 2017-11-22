/*
 * VoronoiPrint.h
 *
 *  Created on: Oct 13, 2011
 *      Author: marchi
 */

#ifndef VORONOIPRINT_H_
#define VORONOIPRINT_H_

class VoronoiPrint {
public:
	static int pGroup;
	static bool bPrintVols;
	static bool bPrintAreas;
	static bool bPrintShell;
	VoronoiPrint();
	virtual ~VoronoiPrint();
};

#endif /* VORONOIPRINT_H_ */
