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
	void WriteIt(std::ofstream &);
	void __extraInit(Topol &,bool);
	void __compShell();
	void __searchNeighs(int ,int);
	void __computeAggregate();
public:

	VoronoiMicelles(Topol &,bool);

	void getData();
	virtual ~VoronoiMicelles();
};
}
#endif /* VORONOIMICELLES_H_ */
