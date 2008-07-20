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


#include <cmath>
#include <iostream>
#include <fstream>

#include "WaveletDecomp.h"

using namespace std;
using namespace blitz;
using namespace bwave;


int main(int argc, char *argv[], char *envp[]) {
	try {
		Array<unsigned char, 3> img;
		if (argc==2) readPNM(argv[1], img);
		else readPNM("testimage.pnm", img);

		typedef int numtype;

		Array<numtype, 3> intimg(img.shape());
		intimg = numtype(1)*img;

		Array<numtype, 3> origData(intimg.shape()); origData=intimg;

		double timer = runTime();
		WaveletDecomp<3> decomp(WL_CDF_2_2, NONSTD_DECOMP, 0);
		decomp.dimSelection(TinyVector<bool, 3>(true, true, false));


		// Speed test:
		for(int i=0; i<100; ++i) {
			decomp.apply(intimg);
			decomp.applyInv(intimg);
		}

		timer = runTime() - timer;
		cout << "Time: " << (double) timer << " s" << endl;
		cout << "MaxDiff (must be << 1): " << max(abs(origData-intimg)) << endl;


		// Coefficient sorting test:

		TinyVector<int, 3> decDepth = decomp.apply(intimg);
		Array<numtype, 3> sepimg(intimg.shape());
		WaveletDecomp<3> decomp2(decomp); decomp2.storageMode(SEPARATED_COEFFS);
		Array<TinyVector<int,3>, 1> idx( decomp2.indices(intimg) );
		for (int i=0; i<idx.rows(); ++i)
			decomp2.coeffs(sepimg, idx(i)) = decomp.coeffs(intimg, idx(i));
		intimg=0;
		for (int i=0; i<idx.rows(); ++i)
			decomp.coeffs(intimg, idx(i)) = decomp2.coeffs(sepimg, idx(i));
		TinyVector<int, 3> recDepth = decomp.applyInv(intimg);
		cout << "MaxDiff (must be << 1): " << max(abs(origData-intimg)) << endl;

		cout << "Decomposition depth = " << decDepth << endl;
		cout << "Recomposition depth = " << recDepth << endl;
		/*cout << "Decomposition indices:" << endl;
		for (int i=0; i<idx.rows(); ++i) {
			cout << "\t" << idx(i) << endl;
		}*/
		

		// Quantisation test:
		decomp.apply(intimg);

		for (int ix=0; ix<idx.rows(); ++ix) {
			TinyVector<int,3> index = idx(ix);
			Array<numtype, 3> coeffs = decomp.coeffs(intimg, index);

			//  obsolete:
			// int level=0, type=0;
			// for (int i=0; i<3; ++i) {
			// 	level = max(level, abs(index(i)));
			// 	if (index(i)<0) ++type;
			// }
			//
			// double thresh = 10;
			// thresh *= pow( sqrt(2.0) , double(type));
			// thresh /= pow( sqrt(2.0) , 2*double(level));

			double thresh = 9 / decomp.normFactor(index);
			numtype ithresh = numtype(thresh);

			for (int i=0; i<coeffs.rows(); ++i)
				for (int j=0; j<coeffs.cols(); ++j)
					for (int k=0; k<coeffs.depth(); ++k)
						coeffs(i,j,k) = abs(coeffs(i,j,k)) > ithresh ? coeffs(i,j,k) : 0;
		}
		

		for (int ix=0; ix<idx.rows(); ++ix)
			decomp2.coeffs(sepimg, idx(ix)) = decomp.coeffs(intimg, idx(ix));

		decomp.applyInv(intimg);
		for (int i=0; i<intimg.rows(); ++i)
			for (int j=0; j<intimg.cols(); ++j)
				for (int k=0; k<intimg.depth(); ++k)
					intimg(i,j,k) = max(numtype(0), min(numtype(255), intimg(i,j,k)) );
					
		cout << "MaxDiff (should be small): " << max(abs(origData-intimg)) << endl;
		cout << "StdErr  (should be small): " << sqrt( 1. * sum((origData-intimg)*(origData-intimg)) / (origData.size()-1) ) << endl;


		// Output image generation:

		Array<unsigned char, 3> outimg(intimg.shape());

		sepimg = sepimg  / 2 + 128;
		outimg = (unsigned char) 1*sepimg;
		writePNM("decomp.pnm", outimg);

		//intimg = abs(origData - intimg);
		outimg = (unsigned char) 1*intimg;
		writePNM("recomp.pnm", outimg);


		// Wavelet form generation:

		WaveletDecomp<1> wldec(WL_CDF_2_2, NONSTD_DECOMP, 0, NESTED_COEFFS);
		Array<double, 1> rwavelet(256); rwavelet=0;
		wldec.coeffs(rwavelet, shape(-6))(1) = 1 / wldec.normFactor(shape(-6));
		wldec.applyInv(rwavelet);
		
		Array<double, 1> awavelet(256); awavelet=0;
		for (int i=0; i<255; ++i) {
			Array<double, 1> wltmp(256); wltmp=0; wltmp(i)=1;
			wldec.apply(wltmp);
			awavelet(i) = wldec.coeffs(wltmp, shape(-6))(1) * wldec.normFactor(shape(-6));
		}

		cout << "sum(a-wavelet * a-wavelet): " << sum(awavelet*awavelet) << endl;
		cout << "sum(r-wavelet * r-wavelet): " << sum(rwavelet*rwavelet) << endl;
		cout << "sum(a-wavelet * r-wavelet): " << sum(awavelet*rwavelet) << endl;

		ofstream wlFormOut("waveletForm.txt");
		for (int i=0; i<awavelet.rows(); ++i) {
			wlFormOut << i << "\t" << awavelet(i) << "\t" << rwavelet(i) << endl;
		}

	}
	catch(std::exception &e) {
		cerr << endl << endl << "Exception: " << e.what() << endl << endl;
	}

	return 0;
}
