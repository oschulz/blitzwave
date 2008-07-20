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


#include "Wavelet.h"
#include "WaveletDecomp.h"

using namespace std;

namespace bwave {


Wavelet::LiftingStep::LiftingStep (Type t, int o, double d, double c0)
	: m_type(t), m_origin(o), m_coeffs(1), m_divisor(d) { m_coeffs[0]=c0; }

Wavelet::LiftingStep::LiftingStep (Type t, int o, double d, double c0, double c1)
	: m_type(t), m_origin(o), m_coeffs(2), m_divisor(d) { m_coeffs[0]=c0; m_coeffs[1]=c1; }

Wavelet::LiftingStep::LiftingStep (Type t, int o, double d, double c0, double c1, double c2)
	: m_type(t), m_origin(o), m_coeffs(3), m_divisor(d) { m_coeffs[0]=c0; m_coeffs[1]=c1; m_coeffs[2]=c2; }

Wavelet::LiftingStep::LiftingStep (Type t, int o, double d, double c0, double c1, double c2, double c3)
	: m_type(t), m_origin(o), m_coeffs(4), m_divisor(d) { m_coeffs[0]=c0; m_coeffs[1]=c1; m_coeffs[2]=c2; m_coeffs[3]=c3; }

Wavelet::LiftingStep::LiftingStep (Type t, int o, double d, double c0, double c1, double c2, double c3, double c4)
	: m_type(t), m_origin(o), m_coeffs(5), m_divisor(d) { m_coeffs[0]=c0; m_coeffs[1]=c1; m_coeffs[2]=c2; m_coeffs[3]=c3; m_coeffs[4]=c4; }

Wavelet::LiftingStep::LiftingStep (Type t, int o, double d, double c0, double c1, double c2, double c3, double c4, double c5)
	: m_type(t), m_origin(o), m_coeffs(6), m_divisor(d) { m_coeffs[0]=c0; m_coeffs[1]=c1; m_coeffs[2]=c2; m_coeffs[3]=c3; m_coeffs[4]=c4; m_coeffs[5]=c5; }


bool Wavelet::LiftingStep::operator==(const LiftingStep &s) const {
	if (!( (type()==s.type()) && equals(origin(), s.origin()) && equals(size(), s.size())
	       && equals(divisor(), s.divisor()) )) return false;
	for (int i=0; i<size(); ++i) if ( !equals(coeff(i), s.coeff(i)) ) return false;
	return true;
}


std::ostream & operator<<(std::ostream &os, const Wavelet::LiftingStep &s) {
	if ( (s.type()==Wavelet::LiftingStep::PRIMAL) || 
		 (s.type()==Wavelet::LiftingStep::DUAL) )
	{
		string src = ( s.type()==Wavelet::LiftingStep::PRIMAL ? "d" : "s" );
		string trg = ( s.type()==Wavelet::LiftingStep::PRIMAL ? "s" : "d" );
		os << ( s.type()==Wavelet::LiftingStep::PRIMAL ? "primal lifting: " : "dual lifting:   " )
		   << trg << "(i) = " << trg << "(i) +";
		if (!Wavelet::equals(s.divisor(),1.0)) os << " 1/" << s.divisor();
		os << " (" ;
		for (int i=0; i<s.size(); ++i) {
			os << (s.coeff(i)<0 ? " - " : (i>0 ? " + " : " ") );
			if (!Wavelet::equals(fabs(s.coeff(i)),1.0)) os << fabs(s.coeff(i)) << " ";
			os << src << "(i" << (s.origin()+i > 0 ? "+" : "");
			if (s.origin()+i != 0) os << s.origin()+i;
			os << ")";
		}
		os << " )";
	}
	else if ( (s.type()==Wavelet::LiftingStep::SSCALE) || 
		 (s.type()==Wavelet::LiftingStep::DSCALE) )
	{
		string trg = ( s.type()==Wavelet::LiftingStep::SSCALE ? "s" : "d" );
		os << trg << "-scaling:      ";
		os << trg << "(i) = " << s.coeff(0);
		if (s.divisor()!=1.) os << "/" << s.divisor();
		os << " * " << trg << "(i)";
	}
	else assert(false);

	return os;
}


bool Wavelet::equals(double x, double y)
	{ return fabs(x-y) < std::numeric_limits<double>::epsilon()*2; }


bool Wavelet::operator==(const Wavelet &w) const {
	if (!( equals(steps(), w.steps()) && equals(normS(), w.normS())
	       && equals(normD(), w.normD()) )) return false;
	for (int i=0; i<steps(); ++i) if (! (liftingStep(i) == w.liftingStep(i) )) return false;
	return true;
}


Wavelet::Wavelet(const std::string &spec, double ns, double nd, const LiftingStep &s0, const LiftingStep &s1)
	: m_spec(spec), m_nS(ns), m_nD(nd), m_steps(2) { m_steps[0]=s0; m_steps[1]=s1; }

Wavelet::Wavelet(const std::string &spec, double ns, double nd, const LiftingStep &s0, const LiftingStep &s1, const LiftingStep &s2)
	: m_spec(spec), m_nS(ns), m_nD(nd), m_steps(3) { m_steps[0]=s0; m_steps[1]=s1; m_steps[2]=s2; }

Wavelet::Wavelet(const std::string &spec, double ns, double nd, const LiftingStep &s0, const LiftingStep &s1, const LiftingStep &s2, const LiftingStep &s3)
	: m_spec(spec), m_nS(ns), m_nD(nd), m_steps(4) { m_steps[0]=s0; m_steps[1]=s1; m_steps[2]=s2; m_steps[3]=s3; }



blitz::Array<double, 1> Wavelet::forwardFkt(unsigned int size,
	int scale, unsigned int trans) const
{
		using namespace blitz;
		WaveletDecomp<1> wldec(*this, NONSTD_DECOMP, ::abs(scale), NESTED_COEFFS);
		Array<double, 1> data(size); data=0;
		for (unsigned int i=0; i+1<size; ++i) {
			Array<double, 1> tmp(size); tmp=0; tmp(i)=1;
			wldec.apply(tmp);
			data(i) = wldec.coeffs(tmp, shape(scale))(trans) * wldec.normFactor(shape(scale));
		}
		return data;
}


blitz::Array<double, 1> Wavelet::inverseFkt(unsigned int size,
	int scale, unsigned int trans) const
{
		using namespace blitz;
		WaveletDecomp<1> wldec(*this, NONSTD_DECOMP, ::abs(scale), NESTED_COEFFS);
		Array<double, 1> data(size); data=0;
		wldec.coeffs(data, shape(scale))(trans) = 1 / wldec.normFactor(shape(scale));
		wldec.applyInv(data);
		return data;
}


	
std::ostream & operator<<(std::ostream &os, const Wavelet &w) {
	os << "Wavelet: " << w.name() << endl;
	// os << "\tCoefficient initialization:" << endl;
	// os << "\t\tLF coefficients s(i) = data(2i)" << endl;
	// os << "\t\tHF coefficients d(i) = data(2i+1)." << endl;
	os << "\tLifting steps:" << endl;
	for (int i=0; i<w.steps(); ++i) os << "\t\t" << w.liftingStep(i) << endl;
	os << "\tNormalisation:" << endl;
	os << "\t\ts(i) = s(i) * " << w.normS() << endl;
	os << "\t\td(i) = d(i) * " << w.normD() << endl;
	return os;
}


const Wavelet WL_CDF_1_1( "CDF(1,1)", sqrt(2.0), sqrt(2.0)/2,
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,    0,  1, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL,  0,  2,  1)
);

const Wavelet WL_CDF_2_2( "CDF(2,2)", sqrt(2.0), sqrt(2.0)/2,
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,    0,  2, -1, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL, -1,  4,  1,  1)
);

const Wavelet WL_CDF_3_1( "CDF(3,1)", 3*sqrt(2.)/2, sqrt(2.)/3,
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL, -1,  3, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,    0,  8, -9, -3),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL,  0,  9,  4)
);

const Wavelet WL_CDF_3_3( "CDF(3,3)", 3*sqrt(2.)/2, sqrt(2.)/3,
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL, -1,  3, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,    0,  8, -9, -3),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL, -1,  36,  3, 16, -3)
);

const Wavelet WL_CDF_4_2( "CDF(4,2)", 2*sqrt(2.), sqrt(2.)/4,
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL, -1,  4, -1, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,    0,  1, -1, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL, -1, 16,  3,  3)
);

const Wavelet WL_CDF_97( "CDF(9,7)", 1.149604398, 1/1.149604398,
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,    0,  1/1.586134342,   -1,  -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL, -1,  1/0.05298011854, -1,  -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,    0,  1/0.8829110762,   1,   1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL, -1,  1/0.4435068522,   1,   1)
);

const Wavelet WL_HAAR( "Haar", sqrt(2.0), sqrt(2.0)/2,
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,    0,  1, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL,  0,  2,  1)
);

const Wavelet WL_LEG_5_3( "LeGall(5,3)", sqrt(2.0), sqrt(2.0)/2,
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,    0,  2, -1, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL, -1,  4,  1,  1)
);

const Wavelet WL_CUBIC_SPLINE( "CubicSpline", sqrt(2.), sqrt(2.)/4,
	Wavelet::LiftingStep(Wavelet::LiftingStep::SSCALE,   0,  1, 2),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL,  -1,  2, -1, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,     0,  2, -1, -1),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL,  -1,  8,  3,  3)
);

const Wavelet WL_D_4( "Daubechies(4)", (sqrt(3.)-1)/sqrt(2.), (sqrt(3.)+1)/sqrt(2.),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL,  0,  1, sqrt(3.)),
	Wavelet::LiftingStep(Wavelet::LiftingStep::DUAL,   -1,  4, -(sqrt(3.)-2.), -sqrt(3.)),
	Wavelet::LiftingStep(Wavelet::LiftingStep::PRIMAL,  1,  1, -1)
);


} // namespace bwave
