// Copyright (C) 2003-2008 Oliver Schulz <oliver.schulz@tu-dortmund.de>

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.


#include <iostream>
#include <fstream>
#include <sys/time.h>

#include "arrayTools.h"

using namespace std;
using namespace blitz;

namespace bwave {


double runTime() {
	timeval t;
	gettimeofday(&t, 0);
	return double(t.tv_sec) + double(t.tv_usec)/double(1000000);
}


void readPNM(const std::string &filename, blitz::Array<unsigned char, 3> &target) {
	std::ifstream source(filename.c_str());

	string format; format+=source.get(); format+=source.get();
	
	int w, h, maxval;
	source >> w >> h >> maxval;
	source.get();

	if (maxval!=255) {
		cerr << "Error: Unsupported max value!" << endl;
		exit(1);
	}
	
	if (format=="P5") {
		if ((target.cols()!=w) || (target.rows()!=h) || (target.depth()!=1))
			target.resize(shape(h, w, 1));
	} else if (format=="P6") {
		if ((target.cols()!=w) || (target.rows()!=h) || (target.depth()!=3))
			target.resize(shape(h, w, 3));
	} else {
		cerr << "Error: Unknown PNM Format \"" << format << "\"!" << endl;
		exit(1);
	}
	
	if (!target.isStorageContiguous()) exit(1);

	source.read((char*)target.dataFirst(), target.size()*sizeof(char));
}


void writePNM(const std::string &filename, blitz::Array<unsigned char, 3> &source) {
	std::ofstream target(filename.c_str());

	if (!source.isStorageContiguous()) exit(1);

	if (source.depth() == 1) target << "P5" << char(10);
	if (source.depth() == 3) target << "P6" << char(10);
	target << source.cols() << " " << source.rows() << char(10) << 255 << char(10);
	target.write ((char*)source.dataFirst(), source.size()*sizeof(char));
}


} // namespace bwave
