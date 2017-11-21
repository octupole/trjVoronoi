/*
 * VoronoiPrint.h
 *
 *  Created on: Oct 13, 2011
 *      Author: marchi
 */

#ifndef VORONOIPRINT_H_
#define VORONOIPRINT_H_
#ifndef FALSE
#  define FALSE   0
#endif
#ifndef TRUE
#  define TRUE    1
#endif

class VoronoiPrint {
public:
	static int pGroup;
	static bool bPrintVols;
	static bool bPrintAreas;
	VoronoiPrint();
	virtual ~VoronoiPrint();
};

#endif /* VORONOIPRINT_H_ */
