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


#ifndef __MAGS_FILTER_H__
#define __MAGS_FILTER_H__


#include <seiscomp3/math/filter.h>
#include <seiscomp3/math/fft.h>
#include <seiscomp3/math/filter/butterworth.h>


template <typename T>
class ButterworthBandstop : public Seiscomp::Math::Filtering::IIR::ButterworthHighLowpass<T> {
	public:
		ButterworthBandstop(int order = 3, double fmin = 0.7, double fmax = 2.0, double fsamp=0);

		virtual void apply(int n, T *inout) {
			for ( int i = 0; i < n; ++i ) {
				T filtered_sample = inout[i];
				Seiscomp::Math::Filtering::IIR::ButterworthHighLowpass<T>::apply(1, &filtered_sample);
				inout[i] -= filtered_sample;
			}
		}
};


template <typename T>
class Envelope : public Seiscomp::Math::Filtering::InPlaceFilter<T> {
	public:
		Envelope(double timeSpan /*sec*/, double fsamp = 0.0)
		: _timeSpan(timeSpan), _fsamp(0.0) {
			if ( fsamp )
				setSamplingFrequency(fsamp);
		}


	// InPlaceFilter interface
	public:
		void setSamplingFrequency(double fsamp) {
			if (_fsamp == fsamp)
				return;

			_fsamp = fsamp;
			_sampleCount = (int)(_fsamp * _timeSpan);
			if ( _sampleCount < 1 ) _sampleCount = 1;
			_sum = 0;
			_index = 0;

			_oocount = 2.0/_sampleCount;

			_buffer.resize(_sampleCount);
			std::fill(_buffer.begin(), _buffer.begin() + _sampleCount, (T)0);
		}

		int setParameters(int n, const double *params) {
			if ( n != 1 ) return 1;
			if ( params[0] <= 0 )
				return -1;

			_timeSpan = params[0];
			return n;
		}

		void apply(int n, T *inout) {
			if ( _fsamp == 0.0 )
				throw Seiscomp::Core::GeneralException("Samplerate not initialized");

			for ( int i = 0; i < n; ++i ) {
				T newValue = inout[i]*inout[i];

				T oldValue = _buffer[_index];
				_buffer[_index] = newValue;

				++_index;

				if ( _index >= _sampleCount )
					_index = 0;

				_sum = _sum + newValue - oldValue;
				inout[i] = (T)sqrt(_sum * _oocount);
			}
		}

		Seiscomp::Math::Filtering::InPlaceFilter<T> *clone() const {
			return new Envelope<T>(_timeSpan, _fsamp);
		}

	private:
		double _timeSpan;
		double _fsamp;
		double _oocount;
		int    _sampleCount;
		int    _index;
		T      _sum;
		std::vector<T> _buffer;
};



void bandpass(std::vector<Seiscomp::Math::Complex> &spec,
              double df, int order, double loFreq, double hiFreq);

void bandstop(std::vector<Seiscomp::Math::Complex> &spec,
              double df, int order, double loFreq, double hiFreq);

void hipass(std::vector<Seiscomp::Math::Complex> &spec,
            double df, int order, double loFreq);

void lopass(std::vector<Seiscomp::Math::Complex> &spec,
            double df, int order, double hiFreq);

void gausslopass(std::vector<Seiscomp::Math::Complex> &spec,
                 double df, double hiFreq);

void hilbert(std::vector<Seiscomp::Math::Complex> &spec);

template <typename T>
void hilbert(const T *x, int n, T *y);

template <typename T>
void magnitude(const T *x, int n, T *y);

template <typename T>
void abs(int n, T *y);

template <typename T>
void envelope(const T *x, int n, T *y);


#endif
