// simpleDecompose.cxx

#include "simpleDecompose.h"
#include <WaveletDecomp.h>

void simpleDecompose(numtype *data, int rows, int cols) {
	using namespace blitz;
	using namespace bwave;

	// Create a blitz array on data, assuming C-style memory layout
	// (data is used directly, not copied):
	GeneralArrayStorage<2> storage;
	storage.ordering() = secondDim, firstDim;
	Array<numtype, 2> array(data, shape(rows, cols), neverDeleteData, storage);

	// Create a CDF(2,2) nonstandard wavelet decomposition in 2 dimensions:
	WaveletDecomp<2> decomp(WL_CDF_2_2, NONSTD_DECOMP);
	
	// Apply the decomposition on the data (in situ):
	decomp.apply(array);
}
