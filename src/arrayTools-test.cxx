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
using namespace bwave;


int main(int argc, char *argv[], char *envp[]) {
	try {
		// Array slicing test:
		Array<int, 2> a(4,3);
		a = 11, 12, 13,
		    21, 22, 23,
		    31, 32, 33,
		    41, 42, 43;

		a.reverseSelf(0); a.reverseSelf(1);
		a.reindexSelf(shape(1,1));

		cout << a << endl;
		cout << slice(a, 0, shape(-55,2)) << endl;


		/*Array<int, 1> b = slice(a, 1, shape(2,-1));
		
		cout << b << endl;

		cout << endl << "Array Extension: " << endl;
		Array<int, 1> ext(b.rows()+5+6);

		ext=0;
		symmExt(b, ext, -1);
		cout << ext << endl;
		
		ext=0;
		symmExt(b, ext, 18);
		cout << ext << endl;
		
		ext=0;
		symmExt(b, ext, -18);
		cout << ext << endl;*/

		// Filter test:

		Array<int, 1> intData(256);
		ifstream intDataIn("intdata.txt");
		for (int i=0; i<intData.rows(); ++i) intDataIn >> intData(i);

		Array<int, 1> filtered(intData.shape()); filtered=5;

		double timer = runTime();
		for (int i=0; i<100000; ++i) {
			int filterData[6] = { 5, -39, 162, 162, -39, 5 };
			GenFilter<int, 6> testGenFilter(-3, filterData, 256);
		
			testGenFilter.apply(intData, filtered);
		}
		timer = runTime() - timer;
		cout << "Time: " << timer << " s" << endl;

		ofstream filteredOut("filtered.txt");
		for (int i=0; i<filtered.rows(); ++i) filteredOut << filtered(i) << endl;

		/*!! // Image I/O test:

		Array<unsigned char, 3> image;
		readPNM("testimage.pnm", image);
		writePNM("out.pnm", image); !!*/
		
		// Array fill test:
		
		{
			cout << endl;

			Array<int, 1> src(6); src = 2,4,5,7,8,9;
			Array<int, 1> trg(13); trg = -1;
			src.reindexSelf(shape(3)); trg.reindexSelf(shape(-2));
			
			cout << src << endl;

			ExtensionMode em = CONSTANT_EXT;
			
			fill(trg, src, -100, em); //cout << trg << endl;
			fill(trg, src, -5, em); //cout << trg << endl;
			fill(trg, src, 3, em); cout << trg << endl;
			fill(trg, src, 7, em); //cout << trg << endl;
			fill(trg, src, 100, em); //cout << trg << endl;
		}
	}
	catch(std::exception &e) {
		cerr << endl << endl << "Exception: " << e.what() << endl << endl;
	}

	return 0;
}
