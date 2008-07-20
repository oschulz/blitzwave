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


#ifndef WAVELET_H
#define WAVELET_H

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <iostream>
#include <blitz/array.h>

#include "arrayTools.h"

namespace bwave {


/// @brief General Wavelet class.
///
/// Instances define specific lifting schemes and therefore specific
/// wavelets. Keep your lifting steps integer compatible if you want
/// do operate on both integer an floating point data-types.
/// Some constant instances of this class are available for instant
/// usage, for example @c WL_HAAR and @c WL_D_4.

class Wavelet {
public:
	/// @brief Compare two numbers within numeric limits.
	/// @return @c true if @c x and @c y and (almost) equal, @c false if not.
	static bool equals(double x, double y);

public:
	/// @brief Superclass for LiftingStep and other possible components
	///        / operations defining a wavelet.

	class WaveletComponent {
	public:
		virtual ~WaveletComponent() { }
	};

	/// @brief Instances define a lifting step in a wavelet lifting scheme.
	///
	/// A lifting step consists of a filter and a divisor.
	/// When applied, the lifting step filters one of two
	/// given data channels (s and d), divides the filtered values by
	/// the divisor and then adds/subtracts them to/from the
	/// data in the other channel. The divisor is necessary
	/// for integer lifting, use only integer filter coefficients
	/// and integer divisor for integer compatibility. You can use such
	/// an integer compatible lifting for both integer and floating point
	/// data-types. Divisors which are powers of two will result
	/// in much faster computation.
	///
	/// In addition to classical primal/dual lifting steps you
	/// can specify scaling steps to scale the d- or s-channel.
	/// However, scaling steps lead to information loss in at least
	/// one transformation direction when operating on integer data.

	class LiftingStep : public WaveletComponent {
	public:
		/// @brief  Constants to model lifting step types.
		enum Type { PRIMAL, DUAL, SSCALE, DSCALE };

	protected:
		Type m_type;
		int m_origin;
		std::vector<double> m_coeffs;
		double m_divisor;

		template<class tp_Type, int tp_size> inline void applyFS
			(blitz::Array<tp_Type, 1> &s, blitz::Array<tp_Type, 1> &d, bool inverse, ExtensionMode be) const;

	public:
		/// @brief  Get type of lifting.
		/// @return PRIMAL for a primal, DUAL for a dual lifting step,
		///         SSCALE, DSCALE for a scaling step.
		Type type() const { return m_type; }

		/// @brief  Get the origin of the lifting step filter.
		/// @return The relative position of the first filter coefficient
		///         for the convolution.
		int origin() const { return m_origin; }

		/// @brief  Get the lifting step filter size.
		/// @return The number of filter coefficients.
		int size() const { return int(m_coeffs.size()); }

		/// @brief  Get a filter coefficient.
		/// @param  i index of the coefficient (i>=0)
		/// @return The i-th filter coefficient.
		double coeff(int i) const { return m_coeffs[i]; }

		/// @brief  Get the divisor.
		/// @return The lifting step divisor.
		double divisor() const { return m_divisor; }
		
		/// @brief  Apply the lifting step on some data.
		/// @param  s  "Low-frequency-" or "Phi-" data channel.
		/// @param  d  "High-frequency-" or "Psi-" data channel.
		/// @param  inverse Do inverse lifting.
		/// @param  be Use this extension mode for boundary handling.
		template<class tp_Type> inline void apply
			(blitz::Array<tp_Type, 1> &s, blitz::Array<tp_Type, 1> &d, bool inverse, ExtensionMode be=CONSTANT_EXT) const;

		/// @brief  Default constructor.
		LiftingStep () {}

		/// @brief  Construct a LiftingStep from it's specification.
		/// @param  t  PRIMAL for a primal, DUAL for a dual lifting step,
		///            SSCALE, DSCALE for a scaling step.
		/// @param  o  Filter origin.
		/// @param  d  Divisor.
		/// @param  c0 First filter coefficient, use @c c1, @c c2, etc. for more.
		///
		/// Use only this constructor (i.e. use only one filter coefficient)
		/// to specify a scaling step.
		LiftingStep (Type t, int o, double d, double c0);
		/// @see    LiftingStep (Type t, int o, double d, double c0)
		LiftingStep (Type t, int o, double d, double c0, double c1);
		/// @see    LiftingStep (Type t, int o, double d, double c0)
		LiftingStep (Type t, int o, double d, double c0, double c1, double c2);
		/// @see    LiftingStep (Type t, int o, double d, double c0)
		LiftingStep (Type t, int o, double d, double c0, double c1, double c2, double c3);
		/// @see    LiftingStep (Type t, int o, double d, double c0)
		LiftingStep (Type t, int o, double d, double c0, double c1, double c2, double c3, double c4);
		/// @see    LiftingStep (Type t, int o, double d, double c0)
		LiftingStep (Type t, int o, double d, double c0, double c1, double c2, double c3, double c4, double c5);
		
		/// @brief  Compare two lifting steps.
		/// @param  s target to compare with
		/// @return @c true if step types, filters and divisors are equal,
		///         @c false if not.
		bool operator==(const LiftingStep &s) const;
		
		/// @brief  Destructor.
		virtual ~LiftingStep() {}
	};

protected:
	std::string m_spec;
	double m_nS, m_nD;
	std::vector<LiftingStep> m_steps;

public:
	/// @brief  Get the name of the wavelet.
	const std::string& name() const { return m_spec; }

	/// @brief  Query if wavelet is integer compatible.
	/// @return @c true if wavelet is integer compatible, @c false if not.
	///
	/// Integer compatibility means that wavelet decomposition and
	/// reconstruction do not use floating point operations on the
	/// data. Wavelet decomposition and reconstruction may still
	/// lead to information loss on integer data if the lifting
	/// scheme contains scaling steps.
	///
	/// Not implemented yet!
	bool integerCompatible() const;

	/// @brief  Query if wavelet and scaling coefficients are normalized after
	///         lifting.
	/// @return @c true if normS() and normD() both equal 1,
	/// @c false if not.
	bool normalized() const;

	/// @brief  Get number of steps in the lifting scheme.
	int steps() const { return int(m_steps.size()); }

	/// @brief  Get a step of the lifting scheme.
	/// @param  i  index of the step to get (0 <= i < steps()).
	const LiftingStep& liftingStep(int i) const { return m_steps[i]; }

	/// @brief  Get the normalization factor for the scaling coefficients.
	///
	/// Wavelet and scaling coefficients are not automatically normalized
	/// after lifting, otherwise integer applications would not be possible.
	double normS() const { return m_nS; }

	/// @brief  Get the normalization factor for the wavelet coefficients.
	/// @see    normS().
	double normD() const { return m_nD; }

	/// @brief  Compare with another wavelet.
	/// @param  w  the wavelet to compare with.
	/// @return @c true if the wavelets are equal, @c false if not.
	bool operator==(const Wavelet &w) const;

	/// @brief  Apply the lifting scheme on the data channels s and d.
	/// @param  s   usually initialised with the even original data values,
	///             contains the scaling coefficients after lifting.
	/// @param  d   usually initialised with the odd original data values,
	///             contains the wavelet coefficients after lifting.
	/// @param  be  extension method to use for boundary handling.
	template<class tp_Type> void lift( blitz::Array<tp_Type, 1> &s,
		blitz::Array<tp_Type, 1> &d, ExtensionMode be=CONSTANT_EXT) const;

	/// @brief  Apply the inverse lifting scheme on the data channels s and d.
	/// @param  s  usually initialised with the scaling coefficients.
	/// @param  d  usually initialised with the wavelet coefficients.
	/// @param  be  extension method to use for boundary handling.
	template<class tp_Type> void invLift( blitz::Array<tp_Type, 1> &s,
		blitz::Array<tp_Type, 1> &d, ExtensionMode be=CONSTANT_EXT) const;

		
	/// @brief  Get sampled values (i.e. the form) of the forward transformation
	/// scaling or wavelet function.
	/// @param  size  number of samples.
	/// @param  scale the scale of the function, use positive values to generate
	///               scaling, negative values to generate wavelet functions.
	/// @param  trans the translation of the function.
	/// @return Array of function values.
	///
	/// Useful for getting the data to plot wavelet and scaling functions. Use
	/// a higher number of samples for better resolution.
	blitz::Array<double, 1> forwardFkt(unsigned int size = 1<<9, int scale = 6,
		unsigned int trans = 3) const;
								  
	/// @brief  Get sampled values (i.e. the form) of the inverse transformation
	/// scaling or wavelet function.
	/// @see forwardFkt(unsigned int size = 1<<9, int scale = 6, unsigned int trans = 3)
	blitz::Array<double, 1> inverseFkt(unsigned int size = 1<<9, int scale = 6,
		unsigned int trans = 3) const;
								  
	/// @brief  Construct an empty wavelet.
	Wavelet() {}

	/// @brief  Construct a wavelet and it's lifting scheme.
	/// @param  name  the name of the wavelet.
	/// @param  ns  normalization factor for scaling coefficients.
	/// @param  nd  normalization factor for wavelet coefficients.
	/// @param  s0  first lifting step.
	/// @param  s1  second lifting step, more steps can be specified.
	///
	/// The normalisation factors can later be queried using normS()
	/// and normD().
	Wavelet(const std::string &name, double ns, double nd, const LiftingStep &s0, const LiftingStep &s1);
	/// @see    Wavelet(const std::string &name, double ns, double nd, const LiftingStep &s0, const LiftingStep &s1)
	Wavelet(const std::string &name, double ns, double nd, const LiftingStep &s0, const LiftingStep &s1, const LiftingStep &s2);
	/// @see    Wavelet(const std::string &name, double ns, double nd, const LiftingStep &s0, const LiftingStep &s1)
	Wavelet(const std::string &name, double ns, double nd, const LiftingStep &s0, const LiftingStep &s1, const LiftingStep &s2, const LiftingStep &s3);

	virtual ~Wavelet() {}
};


/// @brief  Print a description a lifting step.
/// @param  os  the stream to print to.
/// @param  s   the lifting step to print.
std::ostream & operator<<(std::ostream &os, const Wavelet::LiftingStep &s);

/// @brief  Print a description a wavelet.
/// @param  os  the stream to print to.
/// @param  w   the wavelet to print.
std::ostream & operator<<(std::ostream &os, const Wavelet &w);


/// @brief  The Cohen-Daubechies-Feauveau (CDF) (1,1) wavelet,
///         identical to the Haar wavelet.
extern const Wavelet WL_CDF_1_1;

/// @brief  The CDF(2,2) biorthogonal wavelet.
///
/// This wavelet (also known as LeGall(5,3)) is used in
/// the JPEG2000 standard.
extern const Wavelet WL_CDF_2_2;

/// @brief  The CDF(3,1) biorthogonal wavelet.
extern const Wavelet WL_CDF_3_1;

/// @brief  The CDF(3,3) biorthogonal wavelet.
extern const Wavelet WL_CDF_3_3;

/// @brief  The CDF(4,2) biorthogonal wavelet.
extern const Wavelet WL_CDF_4_2;

/// @brief  The CDF(9,7) biorthogonal wavelet (used in JPEG 2000).
///
/// This wavelet (often called Daubechies(9,7) or just 9/7-filter
/// in the literature) is used in the JPEG2000 standard.
extern const Wavelet WL_CDF_97;

/// @brief  The Haar wavelet.
extern const Wavelet WL_HAAR;

/// @brief  The LeGall(5,3) biorthogonal wavelet, identical to CDF(2,2).
///
/// This wavelet is used in the JPEG2000 standard.
extern const Wavelet WL_LEG_5_3;

/// @brief  A cubic-spline biorthogonal wavelet, equals CDF(4,2) except
///         for normalization.
///
/// Lifting includes a scaling step leading to information loss for
///	integer values during recomposition.
extern const Wavelet WL_CUBIC_SPLINE;

/// @brief  The Daubechies D4 wavelet.
extern const Wavelet WL_D_4;


template<class tp_Type, int tp_size> void Wavelet::LiftingStep::applyFS
	(blitz::Array<tp_Type, 1> &s, blitz::Array<tp_Type, 1> &d, bool inverse, ExtensionMode be) const
{
	assert(size() == tp_size);
	tp_Type coeffs[tp_size];
	for (int i=0; i<size(); ++i) coeffs[i] = (tp_Type) coeff(i);
	GenFilter<tp_Type, tp_size> liftFilter(origin(), coeffs, (tp_Type) divisor());
		
	if (type() == PRIMAL) {
		if (!inverse) liftFilter.applyAdd(d, s, be);
		else liftFilter.applySub(d, s, be);
	} else if (type() == DUAL) {
		if (!inverse) liftFilter.applyAdd(s, d, be);
		else liftFilter.applySub(s, d, be);
	} else if (type() == SSCALE) {
		if (!inverse) {
			s *= tp_Type(coeff(0));
			if (divisor()!=1.) s /= tp_Type(divisor()); //!! Add speed-up for divisor()==2^n
		}
		else {
			if (divisor()!=1.) s *= tp_Type(divisor()); //!! Add speed-up for divisor()==2^n
			s /= tp_Type(coeff(0));
		}
	} else if (type() == DSCALE) {
		if (!inverse) {
			d *= tp_Type(coeff(0));
			if (divisor()!=1.) d /= tp_Type(divisor()); //!! Add speed-up for divisor()==2^n
		}
		else {
			if (divisor()!=1.) d *= tp_Type(divisor()); //!! Add speed-up for divisor()==2^n
			d /= tp_Type(coeff(0));
		}
	} else {
		assert(false);
	}
}


template<class tp_Type> void Wavelet::LiftingStep::apply
	(blitz::Array<tp_Type, 1> &s, blitz::Array<tp_Type, 1> &d, bool inverse, ExtensionMode be) const
{
	switch(size()) {
		case 1: applyFS<tp_Type, 1>(s, d, inverse, be); break;
		case 2: applyFS<tp_Type, 2>(s, d, inverse, be); break;
		case 3: applyFS<tp_Type, 3>(s, d, inverse, be); break;
		//case 4: applyFS<tp_Type, 4>(s, d, inverse, be); break;
		//case 5: applyFS<tp_Type, 5>(s, d, inverse, be); break;
		//case 6: applyFS<tp_Type, 6>(s, d, inverse, be); break;
		default: assert(false);
	}
}


template<class tp_Type> void Wavelet::lift( blitz::Array<tp_Type, 1> &s,
	blitz::Array<tp_Type, 1> &d, ExtensionMode be) const
{
	for (int step=0; step<steps(); ++step)
		liftingStep(step).apply(s, d, false, be);
}


template<class tp_Type> void Wavelet::invLift( blitz::Array<tp_Type, 1> &s,
	blitz::Array<tp_Type, 1> &d, ExtensionMode be) const
{
	for (int step=steps()-1; step>=0; --step)
		liftingStep(step).apply(s, d, true, be);
}


} // namespace bwave

#endif
