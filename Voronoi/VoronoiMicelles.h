/*
 * VoronoiMicelles.h
 *
 *  Created on: Nov 19, 2017
 *      Author: marchi
 */

#ifndef VORONOIMICELLES_H_
#define VORONOIMICELLES_H_

#include <Voronoi.h>
namespace Voro{
class VoronoiMicelles: public Voro::Voronoi {
protected:
	void WriteIt(std::ofstream &);
	void __extraInit(Topol &,bool);
	void __compShell();
	void __searchNeighs(int ,int);
	void __computeAggregate();
	VoronoiMicelles(){};
public:

	VoronoiMicelles(ifstream & f){};
	VoronoiMicelles(Topol &,bool);

	void getData();
	virtual ~VoronoiMicelles();
};
}
#endif /* VORONOIMICELLES_H_ */
