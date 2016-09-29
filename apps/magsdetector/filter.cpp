/***************************************************************************
 *   You can redistribute and/or modify this program under the             *
 *   terms of the SeisComP Public License.                                 *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   SeisComP Public License for more details.                             *
 *                                                                         *
 *   Developed by gempa GmbH                                               *
 ***************************************************************************/


#include "filter.h"
#include <math.h>


template <typename T>
ButterworthBandstop<T>::ButterworthBandstop(int order, double fmin, double fmax, double fsamp)
: Seiscomp::Math::Filtering::IIR::ButterworthHighLowpass<T>(order, fmin, fmax, fsamp) {}


void bandpass(std::vector<Seiscomp::Math::Complex> &spec,
              double df, int order, double loFreq, double hiFreq) {
	double freq = df;
	double exp = order*2.0;

	double eflh = pow(loFreq / hiFreq, exp);
	for ( size_t i = 0; i < spec.size(); ++i ) {
		double v = 1.0 + pow(freq / hiFreq, exp) + pow(loFreq / freq, exp) + eflh;
		spec[i] *= 1.0 / sqrt(v);
		freq += df;
	}
}

void bandstop(std::vector<Seiscomp::Math::Complex> &spec,
              double df, int order, double loFreq, double hiFreq) {
	double freq = df;
	double exp = order*2.0;

	double eflh = pow(loFreq / hiFreq, exp);
	for ( size_t i = 0; i < spec.size(); ++i ) {
		double v = 1.0 + pow(freq / hiFreq, exp) + pow(loFreq / freq, exp) + eflh;
		spec[i] *= 1.0 - 1.0 / sqrt(v);
		freq += df;
	}
}

void hipass(std::vector<Seiscomp::Math::Complex> &spec,
            double df, int order, double loFreq) {
	double freq = df;
	double exp = order*2.0;

	for ( size_t i = 0; i < spec.size(); ++i ) {
		double v = 1.0 + pow(loFreq / freq, exp);
		spec[i] *= 1.0 / sqrt(v);
		freq += df;
	}
}

void lopass(std::vector<Seiscomp::Math::Complex> &spec,
            double df, int order, double hiFreq) {
	double freq = df;
	double exp = order*2.0;

	for ( size_t i = 0; i < spec.size(); ++i ) {
		double v = 1.0 + pow(freq / hiFreq, exp);
		spec[i] *= 1.0 / sqrt(v);
		freq += df;
	}
}


void gausslopass(std::vector<Seiscomp::Math::Complex> &spec,
                 double df, double hiFreq) {
	double freq = df;
	for ( size_t i = 0; i < spec.size(); ++i ) {
		double c = (freq/(2*hiFreq));
		spec[i] *= exp(-(c*c));
		freq += df;
	}
}


void hilbert(std::vector<Seiscomp::Math::Complex> &spec) {
	for ( size_t i = 1; i < spec.size()-1; ++i )
		//spec[i] *= Math::Complex(0,1);
		spec[i] = Seiscomp::Math::Complex(-spec[i].imag(), spec[i].real());
}


template <typename T>
void hilbert(const T *x, int n, T *y) {
	std::vector<Seiscomp::Math::Complex> spec;
	Seiscomp::Math::fft(spec,n,x);
	hilbert(spec);
	Seiscomp::Math::ifft(n,y,spec);
}


template <typename T>
void magnitude(const T *x, int n, T *y) {
	int i = n;
	while ( i ) {
		*y = sqrt(*x**x + *y**y);
		++x; ++y; --i;
	}
}


template <typename T>
void abs(int n, T *y) {
	for ( int i = 0; i < n; ++i )
		y[i] = fabs(y[i]);
}


template <typename T>
void envelope(const T *x, int n, T *y) {
	hilbert(x, n, y);
	magnitude(x, n, y);
}


template void hilbert<float>(const float *x, int n, float *y);
template void hilbert<double>(const double *x, int n, double *y);

template void magnitude<float>(const float *x, int n, float *y);
template void magnitude<double>(const double *x, int n, double *y);

template void abs<float>(int n, float *x);
template void abs<double>(int n, double *x);

template void envelope<float>(const float *x, int n, float *y);
template void envelope<double>(const double *x, int n, double *y);

template class ButterworthBandstop<float>;
template class ButterworthBandstop<double>;
