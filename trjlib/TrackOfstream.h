/*
 * TrackOfstream.h
 *
 *  Created on: Nov 22, 2017
 *      Author: marchi
 */

#ifndef TRJLIB_TRACKOFSTREAM_H_
#define TRJLIB_TRACKOFSTREAM_H_
#include <ostream>
#include <fstream>
#include <string>

using std::string;
using std::ofstream;
class TrackOfstream {
	ofstream * out{nullptr};
	std::ios_base::openmode openMode;
public:
	TrackOfstream(string,std::ios_base::openmode,std::ios_base::openmode);
	std::ios_base::openmode gopenMode(){return openMode;}
	virtual ~TrackOfstream();
};

#endif /* TRJLIB_TRACKOFSTREAM_H_ */
