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


#ifndef ARRAYTOOLS_H
#define ARRAYTOOLS_H

#include <cassert>
#include <iostream>
#include <blitz/array.h>

namespace bwave {


/// @brief Array extension mode.
///
/// Used to specify the array extension method for operations
/// which do boundary handling by extension.
///
/// Possible values are: @c ZERO_EXT, for extension with zeros,
/// @c CONSTANT_EXT for boundary value repetition, @c SYMMETRIC_EXT
/// for symmetric extension, @c SYMMETRIC2_EXT for symmetric extension
/// with boundary value doubling and @c CYCLIC_EXT for cyclic extension
/// (array repetition).

enum ExtensionMode {
	ZERO_EXT=0,
	CONSTANT_EXT=1,
	SYMMETRIC_EXT=2,
	SYMMETRIC2_EXT=3,
	CYCLIC_EXT=2
};


template<class tp_Type, int tp_rank>
	blitz::Array<tp_Type,1> slice(blitz::Array<tp_Type, tp_rank> &array,
	int dimension, blitz::TinyVector<int, tp_rank> position);

template<class tp_Type>
	void fill(blitz::Array<tp_Type,1> &trg, const blitz::Array<tp_Type,1>
		&src, int origin=0, ExtensionMode em=ZERO_EXT);


template<class tp_Type> static tp_Type rollR(tp_Type x, tp_Type i) { assert(false); return 0; }
template<> inline int rollR(int x, int i) { return x >> i; }
template<> inline short rollR(short x, short i) { return x >> i; }
template<> inline char rollR(char x, char i) { return x >> i; }

template<class tp_Type> static tp_Type div2BS(tp_Type divisor) { return -1; }
template<> inline int div2BS(int divisor) { return (divisor>0) && !(divisor & (divisor-1)) ? int( log((double)divisor) / log(2.0) ) : -1; }
template<> inline short div2BS(short divisor) { return (divisor>0) && !(divisor & (divisor-1)) ? short( log((double)divisor) / log(2.0) ) : -1; }
template<> inline char div2BS(char divisor) { return (divisor>0) && !(divisor & (divisor-1)) ? char( log((double)divisor) / log(2.0) ) : -1; }


template<class tp_Type, int tp_size> class GenFilter {
protected:
	tp_Type m_coeffs[tp_size];
	int m_origin;
	tp_Type m_divisor;

	template<class tp_Type2, int tp_size2> static tp_Type2 filterSP
		(const GenFilter<tp_Type2, tp_size2> &coeffs, const blitz::Array<tp_Type2,1> &in, int i)
	{
		tp_Type2 result=0;
		for (int j=0; j<tp_size2; ++j) result += coeffs(j)*in(i+j);
		return result;
	}
	
	template<class tp_Type2> static tp_Type2 filterSP
		(const GenFilter<tp_Type2, 1> &coeffs, const blitz::Array<tp_Type2,1> &in, int i)
	{ return coeffs(0)*in(i+0); }

	template<class tp_Type2> static tp_Type2 filterSP
		(const GenFilter<tp_Type2, 2> &coeffs, const blitz::Array<tp_Type2,1> &in, int i)
	{ return coeffs(0)*in(i+0) + coeffs(1)*in(i+1); }

	template<class tp_Type2> static tp_Type2 filterSP
		(const GenFilter<tp_Type2, 3> &coeffs, const blitz::Array<tp_Type2,1> &in, int i)
	{ return coeffs(0)*in(i+0) + coeffs(1)*in(i+1) + coeffs(2)*in(i+2); }

	template<class tp_Type2> static tp_Type2 filterSP
		(const GenFilter<tp_Type2, 4> &coeffs, const blitz::Array<tp_Type2,1> &in, int i)
	{ return coeffs(0)*in(i+0) + coeffs(1)*in(i+1) + coeffs(2)*in(i+2) + coeffs(3)*in(i+3); }

	template<class tp_Type2> static tp_Type2 filterSP
		(const GenFilter<tp_Type2, 5> &coeffs, const blitz::Array<tp_Type2,1> &in, int i)
	{ return coeffs(0)*in(i+0) + coeffs(1)*in(i+1) + coeffs(2)*in(i+2) + coeffs(3)*in(i+3) + coeffs(4)*in(i+4); }

	template<class tp_Type2> static tp_Type2 filterSP
		(const GenFilter<tp_Type2, 6> &coeffs, const blitz::Array<tp_Type2,1> &in, int i)
	{ return coeffs(0)*in(i+0) + coeffs(1)*in(i+1) + coeffs(2)*in(i+2) + coeffs(3)*in(i+3) + coeffs(4)*in(i+4) + coeffs(5)*in(i+5); }

public:
	static void set(tp_Type &target, const tp_Type &value) { target = value; }
	static void inc(tp_Type &target, const tp_Type &value) { target += value; }
	static void dec(tp_Type &target, const tp_Type &value) { target -= value; }
	
	tp_Type operator()(unsigned i) const { return m_coeffs[i]; }
	tp_Type& operator()(unsigned i) { return m_coeffs[i]; }

	tp_Type divisor() const { return m_divisor; }
	tp_Type divisor(int newDiv) { m_divisor=newDiv; return divisor(); }

	template<void op(tp_Type &target, const tp_Type &value)>
		void apply(const blitz::Array<tp_Type,1> &in, blitz::Array<tp_Type,1> &out, ExtensionMode be) const
	{
		using namespace std;
		assert( in.lbound()(0)==0 );
		assert( out.lbound()(0)==0 );
		
		tp_Type bitShift = div2BS(m_divisor);
		
		const int nI=in.rows();
		const int nO=out.rows();
		const int fb = m_origin;
		const int fe = m_origin+tp_size-1;
		const int mainB = max(0, -fb), mainE = min(nI, nO)+min(0, -fe);

		blitz::Array<tp_Type, 1> dummy(max(1, tp_size-1 + max(mainB, nO-mainE) ));

		if (mainB>0) {
			fill(dummy, in, mainB, be);

			int inI=0;
			if (bitShift>=0) for (int outI=0; outI<mainB; ++outI) {
				op( out(outI), rollR(filterSP(*this, dummy, inI), bitShift) );
				++inI;
			} else for (int outI=0; outI<mainB; ++outI) {
				op( out(outI), filterSP(*this, dummy, inI) / m_divisor );
				++inI;
			}
		}
		
		int inI=mainB+fb;
		if (bitShift>=0) for (int outI=mainB; outI<mainE; ++outI) {
			op( out(outI), rollR(filterSP(*this, in, inI), bitShift) );
			++inI;
		} else for (int outI=mainB; outI<mainE; ++outI) {
			op( out(outI), filterSP(*this, in, inI) / m_divisor );
			++inI;
		}

		if (mainE<nO) {
			fill(dummy, in, tp_size-1-nI, be);

			int inI=0;
			if (bitShift>=0) for (int outI=mainE; outI<nO; ++outI) {
				op( out(outI), rollR(filterSP(*this, dummy, inI), bitShift) );
				++inI;
			} else for (int outI=mainE; outI<nO; ++outI) {
				op( out(outI), filterSP(*this, dummy, inI) / m_divisor );
				++inI;
			}
		}
	}

	void apply(const blitz::Array<tp_Type,1> &in, blitz::Array<tp_Type,1> &out, ExtensionMode be=ZERO_EXT) const
		{ apply<GenFilter<tp_Type, tp_size>::set>(in, out, be); }

	void applyAdd(const blitz::Array<tp_Type,1> &in, blitz::Array<tp_Type,1> &out, ExtensionMode be=ZERO_EXT) const
		{ apply<GenFilter<tp_Type, tp_size>::inc>(in, out, be); }

	void applySub(const blitz::Array<tp_Type,1> &in, blitz::Array<tp_Type,1> &out, ExtensionMode be=ZERO_EXT) const
		{ apply<GenFilter<tp_Type, tp_size>::dec>(in, out, be); }

	GenFilter(int origin) : m_origin(origin) {}

	GenFilter(int origin, const tp_Type *x, tp_Type divis)
		: m_origin(origin), m_divisor(divis)
	{
		for (int i=0; i<tp_size; ++i) m_coeffs[i] = x[i];
	}
	
	~GenFilter() {}
};


double runTime();

void readPNM(const std::string &filename, blitz::Array<unsigned char, 3> &target);

void writePNM(const std::string &filename, blitz::Array<unsigned char, 3> &source);



// -- template implementations: ---------------------------------------

template<class tp_Type, int tp_rank>
	blitz::Array<tp_Type,1> slice(blitz::Array<tp_Type, tp_rank> &array,
	int dimension, blitz::TinyVector<int, tp_rank> position)
{
	using namespace blitz;

	GeneralArrayStorage<1> newStorage;
	newStorage.ordering()(0) = firstDim;
	newStorage.ascendingFlag()(0) = array.stride(dimension)>=0;
	newStorage.base()(0) = array.base(dimension);

	position(dimension) = newStorage.ascendingFlag()(0) ?
		array.lbound(dimension) : array.ubound(dimension);
	tp_Type *data = &(array(position));

	TinyVector<int, 1> newShape; newShape(0) = array.shape()(dimension);
	TinyVector<int, 1> newStride; newStride(0) = array.stride(dimension);

	return Array<tp_Type, 1>(data, newShape, newStride, neverDeleteData, newStorage);
}


template<class tp_Type>
	void fill(blitz::Array<tp_Type,1> &trg, const blitz::Array<tp_Type,1>
		&src, int origin, ExtensionMode em)
{
	using namespace blitz;

	// Symmetric and cyclic extension not implemented yet
	assert ( em!= SYMMETRIC_EXT );
	assert ( em!= SYMMETRIC2_EXT );
	assert ( em!= CYCLIC_EXT );

	//const int nSrc = src.rows();
	//const int nTrg = trg.rows();
	const int lSrc = src.lbound()(0);
	const int lTrg = trg.lbound()(0);
	const int uSrc = src.ubound()(0);
	const int uTrg = trg.ubound()(0);

	const int stop = uTrg+1;
	const int start = lTrg;
	const int sSrc = max(start-origin, 0)+lSrc;
	const int stage1 = min(max(origin, start), stop);
	const int stage2 = min(stage1+max(uSrc-sSrc+1,0), stop);

	tp_Type fillConst = (em==CONSTANT_EXT ? src(lSrc) : 0);
	for (int i=start; i<stage1; ++i) trg(i) = fillConst;

	const int shift = sSrc - stage1;
	for (int i=stage1; i<stage2; ++i) trg(i) = src(i+shift);

	fillConst = (em==CONSTANT_EXT ? src(uSrc) : 0);
	for (int i=stage2; i<stop; ++i) trg(i) = fillConst;
}


} // namespace bwave

#endif
