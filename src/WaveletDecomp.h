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


#ifndef WAVELETDECOMP_H
#define WAVELETDECOMP_H

#include <cassert>
#include <iostream>
#include <blitz/array.h>

#include "Wavelet.h"

namespace bwave {


/// @brief  Decomposition type.
///
/// Possible values for this enum are @c STD_DECOMP for
/// standard decompositions and @c NONSTD_DECOMP for non-standard
/// decompositions.

enum DecompType {
	STD_DECOMP=0,
	NONSTD_DECOMP=1
};


/// @brief  Wavelet/Scaling coefficient storage mode.
///
/// @c NESTED_COEFFS specifies the nested coefficients storage
/// mode native to the in-place lifting scheme operations.
/// Coefficients of different type (scale) are stored alternate
/// in memory.
///
/// @c SEPARATED_COEFFS specifies the separated (or isotropic)
/// storage scheme originally established by Mallat. The
/// coefficients of equal type (scale) are stored in contiguous
/// blocks.
///
/// The class WaveletDecomp provides a method @c coeffs(...) to
/// extract coefficients of equal type, so coefficients can be
/// accessed transparently and independent of the storage scheme.

enum CoeffStorage {
	NESTED_COEFFS=0,
	SEPARATED_COEFFS=1
};


/// @brief  Represents a specific wavelet decomposition
///
/// Instances of WaveletDecomp specify the details of
/// a wavelet decomposition, which can then be applied
/// to data of any type but fixed dimensionality.
/// The data dimensionality is set with the template
/// parameter @c tp_rank.

template<int tp_rank> class WaveletDecomp {
protected:
	Wavelet m_wavelet;
	DecompType m_decomp;
	CoeffStorage m_storageMode;
	ExtensionMode m_extMode;
	int m_maxLevel;
	blitz::TinyVector<bool, tp_rank> m_dimSelect;

	template<class tp_Type> inline void trafoStep
		(blitz::Array<tp_Type,tp_rank> &data, int targetDim, bool inverse) const;
		
	template<class tp_Type> inline blitz::TinyVector<int, tp_rank> waveletDecompose
		(blitz::Array<tp_Type,tp_rank> &data, int maxlevel=0) const;
	
	template<class tp_Type> inline blitz::TinyVector<int, tp_rank> waveletRecompose
		(blitz::Array<tp_Type,tp_rank> &data, int maxlevel=0) const;
	
public:
	/// @brief  Get the wavelet used in this decomposition.
	Wavelet wavelet() const { return m_wavelet; }

	/// @brief  Get the type of this decomposition (i.e. standard
	///         or non-standard).
	DecompType decompType() const { return m_decomp; }

	/// @brief  Get the mode in which the wavelet/scaling coefficients
	///         are stored.
	/// @see CoeffStorage
	CoeffStorage storageMode() const { return m_storageMode; }

	/// @brief  Set the mode in which the wavelet/scaling coefficients
	///         are stored.
	/// @param  newCS the new storage mode
	///
	/// Note that lifting-based operations (i.e. @c apply(...) and
	/// @c applyInv(...) only work in @c NESTED_COEFFS mode.
	/// @see CoeffStorage
	CoeffStorage storageMode(CoeffStorage newCS)
		{ return m_storageMode=newCS; }

	/// @brief  Get the current method used for boundary handling / array
	///        extension.
	/// @see ExtensionMode
	ExtensionMode extensionMode() const { return m_extMode; }

	/// @brief  Set the method to be used for boundary handling / array
	///         extension.
	/// @param  newEM  The new array extension method.
	/// @see ExtensionMode
	ExtensionMode extensionMode(ExtensionMode newEM)
		{ return m_extMode=newEM; }

	/// @brief  Query if decomposition will be applied in array dimension @c dim.
	bool dimSelected(int dim) const { return m_dimSelect(dim); }

	/// @brief  Get the dimensions selected for decomposition.
	/// @see dimSelection(blitz::TinyVector<bool, tp_rank> selection)
	blitz::TinyVector<bool, tp_rank> dimSelection() const { return m_dimSelect; }

	/// @brief  Select the dimensions in which to apply the wavelet decomposition.
	/// @param  selection  A boolean vector, set true for the dimensions in
	///                    which to apply decomposition.
	/// 
	/// The wavelet decomposition can be specified to work in all or only
	/// some dimensions (data directions). All possible combinations are valid.
	blitz::TinyVector<bool, tp_rank> dimSelection(blitz::TinyVector<bool, tp_rank> selection)
		{ return m_dimSelect = selection; }

	/// @brief  Apply this wavelet decomposition to a data array.
	/// @param  data  The data to be decomposed.
	/// @return A vector of the decomposition depths in the
	///         different data dimensions.
	///
	/// Does a wavelet decomposition given on the given data.
	/// The data is wavelet-decomposed in-place, storageMode() must
	/// be @c NESTED_COEFFS. Afterwards, data array contains the
	/// wavelet/scaling coefficients, use the @c coeffs(..) method
	/// for coefficient access.
	/// @see coeffs(blitz::Array<tp_Type,tp_rank> &data, blitz::TinyVector<int, tp_rank> indices)
	template<class tp_Type> inline blitz::TinyVector<int, tp_rank> apply
		(blitz::Array<tp_Type,tp_rank> &data) const
	{
		assert(storageMode()==NESTED_COEFFS);
		return waveletDecompose(data, m_maxLevel);
	}

	/// @brief  Apply the inverse wavelet transformation.
	/// @param  data  The wavelet/scaling coefficients.
	/// @return A vector of the recomposition depths in the
	///         different data dimensions.
	///
	/// Does a wavelet recomposition using the given coefficients.
	/// The recomposition is done in-place, storageMode() must
	/// be @c NESTED_COEFFS. Afterwards, the data array contains the
	/// recomposed data.
	template<class tp_Type> inline blitz::TinyVector<int, tp_rank> applyInv
		(blitz::Array<tp_Type,tp_rank> &data) const
	{
		assert(storageMode()==NESTED_COEFFS);
		return waveletRecompose(data, m_maxLevel);
	}


	/// @brief  Get the indices of the coefficients that this
	///         decomposition will/has produce(d) for the given data.
	/// @param  data  original data or transformed data.
	/// @return An array of coefficient scaling indices which occur
	///         in the decomposed data.
	///
	/// There method is useful if you want to iterate over all
	/// wavelet/scaling coefficients in the decomposed data.
	/// Since wavelet decomposition is done in-place, it does not
	/// matter whether this method is called before or after decomposition
	/// of the data, since the data extents (which determine the possible
	/// indices will stay the same.
	///
	/// @see coeffs(blitz::Array<tp_Type,tp_rank> &data, blitz::TinyVector<int, tp_rank> indices)
	template<class tp_Type>
		blitz::Array< blitz::TinyVector<int, tp_rank>, 1> indices(
			blitz::Array<tp_Type,tp_rank> &data) const;

	/// @brief  Get all wavelet/scaling coefficients of the same type
	///         (i.e. corresponding to basis functions differing only by
	///         translation)
	/// @param  data     The decomposed data to extract the coefficients from.
	/// @param  indices  An index vector specifying the coefficient type to
	///                  extract.
	/// @return  An array of the coefficients. This not a copy of the
	///          relevant parts of @c data, but a reference to them.
	///
	/// Wavelet decomposition of multi-dimensional, non-square data array
	/// produces coefficients coefficients corresponding to lots of
	/// scaled and translated basis functions, which are themselves mixed
	/// products of the basic scaling and wavelet function. Coefficients
	/// which differ only by the translation of their basic function
	/// are said here to be of the same @em type. These coefficient types
	/// are designated here by an index vector. It's components specify the
	/// factory of the basis function in the respective dimensions. A
	/// positive value specifies a scaling function, a negative one a
	/// wavelet function. The component's absolute value specifies the scale
	/// as an exponent of 2.
	///
	/// Example: Let Phi_i be the scaling, Psi_i be the wavelet function
	/// of scale i. Let the data have three dimensions x,y,z.
	/// Then the coefficients corresponding to the basis function
	/// Psi_4(x)*Phi_3(y)*Psi_3(z) would have the index vector (-4, 3, -3).
	///
	/// The bases function Phi_0 can be said to correspond to original,
	/// non-decomposed data. Therefore, if data is not decomposed in
	/// some dimensions, the only valid index vector component in these
	/// dimensions is 0.
	///
	/// @note This method does @em not check if the requested coefficient
	/// type exists for this decomposition and the given coefficient data.
	/// The method @c indices(...) can be used to get a list of valid index
	/// vectors.
	template<class tp_Type>
		blitz::Array<tp_Type,tp_rank> coeffs(blitz::Array<tp_Type,tp_rank> &data,
			blitz::TinyVector<int, tp_rank> indices) const;
	
	/// @brief  Get the normalization factor for the specified coefficient type.
	/// @param  indices  The index vector specifying the coefficient type
	/// @return The factor needed to normalize the coefficients.
	///
	/// The Wavelet and scaling coefficients are not automatically normalized
	/// after decomposition, otherwise integer applications would not be possible.
	/// Also, normalisation is unneeded or even unwanted for some applications.
	double normFactor(blitz::TinyVector<int, tp_rank> indices) const {
		blitz::TinyVector<double, 2> baseNF(wavelet().normS(), wavelet().normD());
		double norm=1;
		for (int i=0; i<tp_rank; ++i) {
			norm *= indices(i)>=0 ? pow(baseNF(0), indices(i))
			                      : pow(baseNF(0), -indices(i)-1) * baseNF(1);
		}
		return norm;
	}


	/// @brief  Construct a new wavelet decomposition.
	/// @param  wavelet  Wavelet to use in the decomposition.
	/// @param  decomp   Decomposition type.
	/// @param  maxlevel Maximum decomposition depth
	/// @param  cs       Coefficient storage mode
	/// @param  em       Internal array extension / boundary handling method
	///
	/// The decomposition can directly be applied to data after construction,
	/// no further preparation is necessary.
	///
	/// @see apply(blitz::Array<tp_Type,tp_rank> &data)
	/// @see applyInv(blitz::Array<tp_Type,tp_rank> &data)
	WaveletDecomp(Wavelet wl, DecompType decomp=NONSTD_DECOMP,
	              int maxlevel=0, CoeffStorage cs=NESTED_COEFFS, ExtensionMode em=CONSTANT_EXT)
		: m_wavelet(wl), m_decomp(decomp), m_storageMode(cs), m_extMode(em),
		  m_maxLevel(maxlevel)
	  { m_dimSelect = true; }
};


template<int tp_rank> template<class tp_Type>
	inline void WaveletDecomp<tp_rank>::trafoStep
	(blitz::Array<tp_Type,tp_rank> &data, int targetDim, bool inverse) const
{
	using namespace blitz;
	assert( targetDim < tp_rank );
	TinyVector<int, tp_rank> stride=1, lboundS=data.lbound(), lboundD=data.lbound();
	stride(targetDim)=2;
	lboundD(targetDim)+=1;

	TinyVector<int, tp_rank> uboundS, uboundD;
	for (int dim=0; dim<tp_rank; ++dim) {
		uboundS(dim) = lboundS(dim) + (data.ubound()(dim)-lboundS(dim)) / stride(dim) * stride(dim);
		uboundD(dim) = lboundD(dim) + (data.ubound()(dim)-lboundD(dim)) / stride(dim) * stride(dim);
	}

	StridedDomain<tp_rank> subsetS(lboundS, uboundS, stride);
	StridedDomain<tp_rank> subsetD(lboundD, uboundD, stride);
	blitz::Array<tp_Type,tp_rank> s( data(subsetS) );
	blitz::Array<tp_Type,tp_rank> d( data(subsetD) );

	int n = data.size() / data.extent(targetDim);

	TinyVector<int, tp_rank> position = data.lbound();
	int lastDim = targetDim!=tp_rank-1 ? tp_rank-1 : tp_rank-2;

	for (int i=0; i<n; ++i) {
		Array<tp_Type, 1> sliceS( slice(s, targetDim, position) );
		Array<tp_Type, 1> sliceD( slice(d, targetDim, position) );
		if (!inverse) m_wavelet.lift(sliceS, sliceD, m_extMode);
		else m_wavelet.invLift(sliceS, sliceD, m_extMode);

		if (i < n-1) {
			position(lastDim) += 1;
			for (int d=lastDim; d>0; --d) {
				if ( (d!=targetDim) && (position(d)>data.ubound(d)) ) {
					position(d)=data.lbound(d);
					++position(targetDim!=d-1 ? d-1 : d-2);
				}
			}
		}
	}

}


template<int tp_rank> template<class tp_Type>
	inline blitz::TinyVector<int, tp_rank>
	WaveletDecomp<tp_rank>::waveletDecompose
	(blitz::Array<tp_Type,tp_rank> &data, int maxlevel) const
{
	using namespace blitz;

	TinyVector<int, tp_rank> depth=0;

	switch (decompType()) {
		case STD_DECOMP: assert(false); break;
		case NONSTD_DECOMP: {
			for (int dim=0; dim<tp_rank; ++dim)
				if ( (data.extent(dim)>1) && dimSelected(dim) ) {
					trafoStep(data, dim, false);
					depth(dim) += 1;
				}
					
			bool descend = false;
			for (int dim=0; dim<tp_rank; ++dim)
				if ((data.extent(dim)>2) && dimSelected(dim)) descend = true;

			if ( descend && ((maxlevel==0) || (maxlevel>1)) ) {
				TinyVector<int, tp_rank> stride, lbound, ubound;
				for (int dim=0; dim<tp_rank; ++dim) {
					stride(dim) = dimSelected(dim) ? 2 : 1;
					lbound(dim) = data.lbound()(dim);
					ubound(dim) = lbound(dim) + (data.ubound()(dim)-lbound(dim)) / stride(dim) * stride(dim);
				}
				StridedDomain<tp_rank> subsetS(lbound, ubound, stride);
				blitz::Array<tp_Type,tp_rank> sub(data(subsetS));
		
				depth += waveletDecompose(sub, maxlevel>0 ? maxlevel-1 : 0);
			}
			
			break;
		}
		default: assert(false);
	}
	
	return depth;
}


template<int tp_rank> template<class tp_Type>
	inline blitz::TinyVector<int, tp_rank>
	WaveletDecomp<tp_rank>::waveletRecompose
	(blitz::Array<tp_Type,tp_rank> &data, int maxlevel) const
{
	using namespace blitz;

	TinyVector<int, tp_rank> depth=0;

	switch (decompType()) {
		case STD_DECOMP: assert(false); break;
		case NONSTD_DECOMP: {
			bool descend = false;
			for (int dim=0; dim<tp_rank; ++dim)
				if ((data.extent(dim)>2) && dimSelected(dim)) descend = true;

			if ( descend && ((maxlevel==0) || (maxlevel>1)) ) {
				TinyVector<int, tp_rank> stride, lbound, ubound;
				for (int dim=0; dim<tp_rank; ++dim) {
					stride(dim) = dimSelected(dim) ? 2 : 1;
					lbound(dim) = data.lbound()(dim);
					ubound(dim) = lbound(dim) + (data.ubound()(dim)-lbound(dim)) / stride(dim) * stride(dim);
				}
				StridedDomain<tp_rank> subsetS(lbound, ubound, stride);
				blitz::Array<tp_Type,tp_rank> sub(data(subsetS));
		
				depth += waveletRecompose(sub, maxlevel>0 ? maxlevel-1 : 0);
			}
		
			for (int dim=tp_rank-1; dim>=0; --dim)
				if ( (data.extent(dim)>1) && dimSelected(dim) ) {
					trafoStep(data, dim, true);
					depth(dim) -= 1;
				}
			
			break;
		}
		default: assert(false);
	}	

	return depth;
}


template<int tp_rank> template<class tp_Type>
	blitz::Array< blitz::TinyVector<int, tp_rank>, 1>
	WaveletDecomp<tp_rank>::indices(
		blitz::Array<tp_Type,tp_rank> &data) const
{
	using namespace blitz;
	TinyVector<int, tp_rank> decompDepth;
	Array< TinyVector<int, tp_rank>, 1> idx;

	for (int dim=0; dim<tp_rank; ++dim) {
		// log(...-0.1) for numerical stability
		int val = (int) ceil(log((double)data.extent(dim)-0.1) / log(2.0));
		decompDepth(dim) = dimSelected(dim) ?
			(m_maxLevel>0 ? min(m_maxLevel, val) : val) : 0;
	}

	switch (decompType()) {
		case STD_DECOMP: assert(false); break;
		case NONSTD_DECOMP: {
			int maxLevel = max(decompDepth);
			idx.resize(max(1, maxLevel*(1<<tp_rank)));
			idx(0) = 0;
			int begin=0, end=0;
			for (int level=1; level<=maxLevel; ++level) {
				for (int dim=0; dim<tp_rank; ++dim)
					if (decompDepth(dim)>=level)
				{
					int offset = end-begin+1;
					for (int i=begin; i<=end; ++i) {
						idx(i+offset) = idx(i);
						idx(i)(dim) = -level;
						idx(i+offset)(dim) = level;
					}
					end += offset;
				}
				begin = end;
			}
			idx.resizeAndPreserve(end+1);
			break;
		}
		default: assert(false);
	}	
	return idx;
}


template<int tp_rank> template<class tp_Type>
	blitz::Array<tp_Type,tp_rank> WaveletDecomp<tp_rank>::coeffs(
		blitz::Array<tp_Type,tp_rank> &data,
		blitz::TinyVector<int, tp_rank> indices) const
{
	using namespace blitz;

	Array<tp_Type,tp_rank> c;
	TinyVector<int, tp_rank> ubound, lbound, stride;
	
	for (int dim=0; dim<tp_rank; ++dim) {
		int i = indices(dim);
		// Remark: i>=0 specifies phi(i) coeff,
		//          i<0 specifies psi(-i) coeff;
		int n = data.extent(dim);
		
		if (storageMode()==NESTED_COEFFS) {
			// Remark: coefficient location in data(l)
			//         for phi(i):       0 <= l <= n-1, stride = 2^i
			//         for psi(i): 2^(i-1) <= l <= n-1, stride = 2^i
			stride(dim) = i>=0 ? 1 << i : 1 << -i;
			lbound(dim) = data.lbound()(dim) + ( i>=0 ? 0 : 1 << (-i-1) );
			ubound(dim) = lbound(dim) + (data.ubound()(dim)-lbound(dim)) / stride(dim) * stride(dim);
		} else if (storageMode()==SEPARATED_COEFFS) {
			// Remark: coefficient location in data(l)
			//         for phi(i):           0 <= l <= ceil(n/2^i)-1
			//         for psi(i): ceil(n/2^i) <= l <= ceil(n/2^(i-1))-1
			// Remark: ceil(n/2^i) = (n + (1<<i) - 1) >> i;
			stride(dim) = 1;
			lbound(dim) = data.lbound()(dim) +  ( i>=0 ? 0 : (n + (1<<-i) - 1) >> -i );
			ubound(dim) = data.lbound()(dim) +  ( i>=0 ? ( (n + (1<<i) - 1) >> i) - 1
													   : ( (n + (1<<(-i-1)) - 1) >> (-i-1) ) - 1 );
		} else assert(false);

		assert( lbound(dim) <= ubound(dim));
		assert( lbound(dim) >= data.lbound(dim));
		assert( ubound(dim) <= data.ubound(dim));
	}

	StridedDomain<tp_rank> subset(lbound, ubound, stride);
	c.reference(data(subset));

	return c;
}


} // namespace bwave

#endif
