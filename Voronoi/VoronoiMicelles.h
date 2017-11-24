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
	virtual void WriteIt(std::ofstream &);
	virtual void __extraInit(Topol &,bool);
	virtual void __compShell();
	virtual void __searchNeighs(int ,int);
public:

	VoronoiMicelles(Topol &,bool);

	void getData();
	virtual ~VoronoiMicelles();
};
}
#endif /* VORONOIMICELLES_H_ */
