/*
 * TrackOfstream.cpp
 *
 *  Created on: Nov 22, 2017
 *      Author: marchi
 */

#include "TrackOfstream.h"

TrackOfstream::TrackOfstream(string filename, std::ios_base::openmode mode, std::ios_base::openmode inOut): openMode{mode}{
	out=new ofstream(filename,mode|inOut);
}

TrackOfstream::~TrackOfstream() {
	if(out) out->close();
}

