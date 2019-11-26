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

#define SEISCOMP_COMPONENT MAGSDETECTOR

#include <seiscomp3/logging/log.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/core/arrayfactory.h>
#include <seiscomp3/core/typedarray.h>
#include <seiscomp3/io/archive/xmlarchive.h>
#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/io/records/mseedrecord.h>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/filter/butterworth.h>

#include <seiscomp3/datamodel/eventparameters.h>
#include <seiscomp3/datamodel/event.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/magnitude.h>
#include <seiscomp3/datamodel/stationmagnitude.h>

#include <seiscomp3/utils/files.h>

#include <iomanip>

#include "detector.h"
#include "filter.h"
#include "crosscorr.h"

#define VERSION "2014.016"
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
using namespace Seiscomp;

#define NEW_OPT(var, ...) addOption(&var, __VA_ARGS__)
#define NEW_OPT_CLI(var, ...) addOption(&var, NULL, __VA_ARGS__)

namespace {


bool checkGain(DataModel::Stream *s) {
	try {
		double gain = s->gain();
		if ( gain == 0.0 ) return false;
	}
	catch ( ... ) {
		return false;
	}

	return true;
}


bool merge(Detector::Trace &trace, RecordSequence *seq) {
	if ( seq->empty() ) {
		SEISCOMP_ERROR("No data");
		return false;
	}

	RecordCPtr first = seq->front();
	RecordCPtr last;
	double samplingFrequency = first->samplingFrequency();
	Core::TimeSpan maxAllowedGap, maxAllowedOverlap;

	maxAllowedGap = Core::TimeSpan((double)(0.5 / samplingFrequency));
	maxAllowedOverlap = Core::TimeSpan((double)(-0.5 / samplingFrequency));

	trace.setNetworkCode(first->networkCode());
	trace.setStationCode(first->stationCode());
	trace.setLocationCode(first->locationCode());
	trace.setChannelCode(first->channelCode());

	trace.setStartTime(first->startTime());
	trace.setSamplingFrequency(samplingFrequency);

	Array::DataType datatype = first->data()->dataType();
	ArrayPtr arr = ArrayFactory::Create(datatype, datatype, 0, NULL);

	RecordSequence::const_iterator it;
	for ( it = seq->begin(); it != seq->end(); ++it ) {
		RecordCPtr rec = *it;
		if ( rec->samplingFrequency() != samplingFrequency ) {
			SEISCOMP_WARNING("%s.%s.%s.%s: record sampling frequencies are not consistent: %f != %f",
			                 trace.networkCode().c_str(),
			                 trace.stationCode().c_str(),
			                 trace.locationCode().c_str(),
			                 trace.channelCode().c_str(),
			                 samplingFrequency, rec->samplingFrequency());
			return false;
		}

		// Check for gaps and overlaps
		if ( last ) {
			Core::TimeSpan diff = rec->startTime()-last->endTime();
			if ( diff > maxAllowedGap ) {
				SEISCOMP_WARNING("%s.%s.%s.%s: gap detected of %d.%06ds",
				                 trace.networkCode().c_str(),
				                 trace.stationCode().c_str(),
				                 trace.locationCode().c_str(),
				                 trace.channelCode().c_str(),
				                 (int)diff.seconds(), (int)diff.microseconds());
				return false;
			}

			if ( diff < maxAllowedOverlap ) {
				SEISCOMP_WARNING("%s.%s.%s.%s: overlap detected of %fs",
				                 trace.networkCode().c_str(),
				                 trace.stationCode().c_str(),
				                 trace.locationCode().c_str(),
				                 trace.channelCode().c_str(),
				                 (double)diff);
				return false;
			}
		}

		arr->append( (Array *)(rec->data()));

		last = rec;
	}

	trace.setData(arr.get());

	return true;
}


bool trim(Detector::Trace &trace, const Core::TimeWindow &tw) {
	int ofs = (int)(double(tw.startTime() - trace.startTime())*trace.samplingFrequency());
	int samples = (int)(tw.length()*trace.samplingFrequency());

	// Not enough data at start of time window
	if ( ofs < 0 ) {
		SEISCOMP_DEBUG("%s: need %d more samples in past",
		               trace.streamID().c_str(), -ofs);
		return false;
	}

	// Not enough data at end of time window
	if ( ofs+samples > trace.data()->size() ) {
		SEISCOMP_DEBUG("%s: need %d more samples past the end",
		               trace.streamID().c_str(), trace.data()->size()-samples-ofs);
		return false;
	}

	ArrayPtr sliced = trace.data()->slice(ofs, ofs+samples);
	/*
	SEISCOMP_DEBUG("%s: cut %d samples at front, cut %d samples at back, size %d, newsize %d",
	               trace.streamID().c_str(),
	               cutFront, cutBack, trace.data()->size(), sliced->size());
	*/

	trace.setStartTime(tw.startTime());
	trace.setData(sliced.get());

	/*
	SEISCOMP_DEBUG("%s: cut to time window %s ~ %s",
	               trace.streamID().c_str(),
	               trace.startTime().iso().c_str(), trace.endTime().iso().c_str());
	*/

	return true;
}


bool resample(Detector::Trace &trace, int sf, bool average) {
	// Resampling is disabled
	if ( sf <= 0 ) return true;
	if ( trace.samplingFrequency() == sf ) return true;

	DoubleArray *data = DoubleArray::Cast(trace.data());
	double step = trace.samplingFrequency() / sf;

	int w = average?step*0.5 + 0.5:0;
	int i = 0;
	double fi = 0.0;
	int cnt = data->size();
	if ( w <= 0 ) {
		while ( fi < cnt ) {
			(*data)[i++] = (*data)[(int)fi];
			fi += step;
		}
	}
	else {
		while ( fi < cnt ) {
			int ci = (int)fi;
			double scale = 1.0;
			double v = (*data)[ci];

			for ( int g = 1; g < w; ++g ) {
				if ( ci >= g ) {
					v += (*data)[ci-g];
					scale += 1.0;
				}

				if ( ci+g < cnt ) {
					v += (*data)[ci+g];
					scale += 1.0;
				}
			}

			v /= scale;

			(*data)[i++] = v;
			fi += step;
		}
	}

	data->resize(i);
	trace.setSamplingFrequency((double)sf);
	trace.dataUpdated();

	return true;
}


void correct(Detector::Trace &trace, double gain) {
	DoubleArray *data = DoubleArray::Cast(trace.data());
	double scale = 1.0 / gain;
	int cnt = data->size();
	for ( int i = 0; i < cnt; ++i )
		(*data)[i] *= scale;
}


template <typename T>
void log(Seiscomp::TypedArray<T> &data) {
	int n = data.size();
	T *d = data.typedData();

	for ( int i = 0; i < n; ++i ) {
		if ( d[i] < 0 )
			d[i] = -::log(-d[i]);
		else if ( d[i] > 0 )
			d[i] = ::log(d[i]);
	}
}



template <typename T>
T norm(int n, const T *d, T &sum, int lag) {
	T norm = 0.0;

	if ( lag < 0 ) {
		n += lag;
		lag = -lag;
	}

	for ( int i = lag; i < n; ++i ) norm += d[i]*d[i];
	sum = norm;
	return sqrt(norm);
}


template <typename T>
void costaper(int n, T *inout, int istart, int iend, int estart, int eend) {
	int taperLength = iend - istart;

	for ( int i = 0; i < istart; ++i )
		inout[i] = 0;

	for ( int i = 0; i < taperLength; ++i ) {
		double frac = double(i)/taperLength;
		inout[istart+i] *= 0.5*(1-cos(M_PI*frac));
	}

	taperLength = eend - estart;

	for ( int i = 0; i < taperLength; ++i ) {
		double frac = double(i)/taperLength;
		inout[estart+i] *= 0.5*(1+cos(M_PI*frac));
	}

	for ( int i = eend; i < n; ++i )
		inout[i] = 0;
}


struct Norm {
	Norm(double scale_syn, double scale_data)
	: scaleSyn(scale_syn), scaleData(scale_data) {}

	double  scaleSyn;
	double  scaleData;
};


struct DotProduct : Norm {
	DotProduct(double scale_syn = 1.0, double scale_data = 1.0)
	: Norm(scale_syn, scale_data), dDOTs(0) {}

	double  dDOTs;

	void update(double syn, double data) {
		double s = syn*scaleSyn;
		double d = data*scaleData;
		dDOTs += d*s;
	}

	DotProduct &operator+=(const DotProduct &other) {
		dDOTs += other.dDOTs;
		return *this;
	}

	DotProduct &scale(double s) {
		dDOTs *= s;
		return *this;
	}

	double operator()() const {
		return dDOTs;
	}
};


typedef DotProduct FitNorm;


DataModel::Origin *findOrigin(DataModel::EventParameters *ep, const std::string &id) {
	for ( size_t i = 0; i < ep->originCount(); ++i )
		if ( ep->origin(i)->publicID() == id )
			return ep->origin(i);
	return NULL;
}


DataModel::Pick *findPick(DataModel::EventParameters *ep, const std::string &id) {
	for ( size_t i = 0; i < ep->pickCount(); ++i )
		if ( ep->pick(i)->publicID() == id )
			return ep->pick(i);
	return NULL;
}


DataModel::Magnitude *findMagnitude(DataModel::Origin *org, const std::string &id) {
	for ( size_t i = 0; i < org->magnitudeCount(); ++i )
		if ( org->magnitude(i)->publicID() == id )
			return org->magnitude(i);
	return NULL;
}


}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::RTTrace::setSamplingFrequency(double sf) {
	_samplingFrequency = sf;
	_maxAllowedGap = Core::TimeSpan((double)(0.5 / _samplingFrequency));
	_maxAllowedOverlap = Core::TimeSpan((double)(-0.5 / _samplingFrequency));
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::RTTrace::setTimeSpan(const Seiscomp::Core::TimeSpan &ts) {
	size_t n = _samplingFrequency * (double)ts + 0.5;
	_data.setCapacity(n);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::RTTrace::reset() {
	_data.clear();
	_startTime = _endTime = Core::Time();
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::RTTrace::feed(Seiscomp::Record *rec) {
	DoubleArrayPtr tmp;
	const DoubleArray *data = DoubleArray::ConstCast(rec->data());
	if ( data == NULL ) {
		tmp = DoubleArray::Cast(rec->data()->copy(Array::DOUBLE));
		data = tmp.get();
	}

	if ( _samplingFrequency != rec->samplingFrequency() ) {
		SEISCOMP_WARNING("%s: inconsistent sampling frequency: %f != %f: ignoring record",
		                 rec->streamID().c_str(), _samplingFrequency,
		                 rec->samplingFrequency());
		return;
	}

	// First record
	if ( !_startTime.valid() ) {
		_startTime = rec->startTime();
		_endTime = rec->endTime();
		_data.clear();
		_data.append(data->typedData(), data->size(), 1.0 / _gain);
		return;
	}

	Core::TimeSpan diff = rec->startTime()-_endTime;
	if ( diff > _maxAllowedGap ) {
		SEISCOMP_WARNING("%s: detected gap of %fs: reset buffer",
		                 rec->streamID().c_str(), (double)diff);
		reset();
		_startTime = rec->startTime();
		_endTime = rec->endTime();
		_data.clear();
		_data.append(data->typedData(), data->size(), 1.0 / _gain);
		return;
	}

	if ( diff < _maxAllowedOverlap ) {
		SEISCOMP_WARNING("%s: detected overlap of %fs: reset buffer",
		                 rec->streamID().c_str(), (double)diff);
		reset();
		_startTime = rec->startTime();
		_endTime = rec->endTime();
		_data.clear();
		_data.append(data->typedData(), data->size(), 1.0 / _gain);
		return;
	}

	_data.append(data->typedData(), data->size(), 1.0 / _gain);

	_endTime = rec->endTime();
	_startTime = _endTime - Core::TimeSpan((double)_data.size() / _samplingFrequency);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::RTTrace::extract(Trace &trace, const Seiscomp::Core::TimeWindow &tw) const {
	trace.setSamplingFrequency(_samplingFrequency);
	DoubleArrayPtr data = new DoubleArray;

	if ( _startTime > tw.startTime() ) return false;
	if ( _endTime < tw.endTime() ) return false;

	trace.setStartTime(tw.startTime());

	size_t samples = (int)(tw.length()*_samplingFrequency);
	size_t ofs = (int)(double(tw.startTime() - _startTime)*_samplingFrequency);

	samples = std::min(samples, _data.size()-ofs);

	/*
	trace.setStartTime(_startTime);
	size_t ofs = 0;
	size_t samples = _data.size();
	*/

	data->resize(samples);
	_data.copy(data->typedData(), ofs, samples);
	//SEISCOMP_DEBUG("copy %d samples", (int)samples);
	//SEISCOMP_DEBUG("st = %s, start = %d, end = %d", _startTime.iso().c_str(), (int)_data._start_index(), (int)_data._end_index());

	trace.setData(data.get());

	return true;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Detector::Config::Config() {
	offline = false;
	test = false;

	processingInterval = 0;
	bufferSize = 600;

	waveParams.filterOrder = 4;
	waveParams.filterFreqs = FilterFreqs(10, 40);
	waveParams.filterBandStop = false;

	waveParams.envelopeEnable = true;
	waveParams.envelopeSamplingFrequency = 20;
	waveParams.envelopeDownsampleAverage = false;
	waveParams.envelopeHiFreq = 0;
	waveParams.envelopeAcausal = false;

	waveParams.processingAcausal = false;
	waveParams.processingLog = false;

	detectorThreshold = 0.55;
	detectorOffThreshold = -1;
	detectorChannelThreshold = 0.55;
	detectorTimeWindow = 2.0;
	detectorMaxTimeLag = 0;
	detectorMinStationRatio = 0;
	detectorMinChannelRatio = 0;
	detectorProcessingTimeWindow = 0;
	detectorPublicationTimeout = 10;
	detectorMaximumStepFrequency = 0;

	dumpFitFunction = 0;
	dumpFitChannels = 0;
	dumpWaveforms = false;
	dumpTriggerFit = false;
	dumpTriggerWaveforms = false;
	maximumLatency = 10.0;

	normalizationString = "total";
	showProgress = false;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::MasterEvent::read(const Seiscomp::Config::Config &cfg,
                                 DataModel::DatabaseQuery *db,
                                 const std::string &prefix) {
	doNotSendOrigins = false;
	cfg.getBool(doNotSendOrigins, prefix + "negative");

	DataModel::EventPtr  evt;
	DataModel::OriginPtr origin;
	DataModel::MagnitudePtr mag;

	try {
		std::string baseEventID = cfg.getString(prefix + "baseID");
		evt = DataModel::Event::Cast(db->getObject(DataModel::Event::TypeInfo(), baseEventID));
		if ( !evt ) {
			SEISCOMP_ERROR("%s.baseID: event not found in database: %s",
			               prefix.c_str(), baseEventID.c_str());
			return false;
		}

		origin = DataModel::Origin::Cast(db->getObject(DataModel::Origin::TypeInfo(), evt->preferredOriginID()));
		if ( !origin ) {
			SEISCOMP_ERROR("%s.baseID: preferred origin not found in database: %s",
			               prefix.c_str(), evt->preferredOriginID().c_str());
			return false;
		}

		mag = DataModel::Magnitude::Cast(db->getObject(DataModel::Magnitude::TypeInfo(), evt->preferredMagnitudeID()));
		if ( !mag ) {
			SEISCOMP_ERROR("%s.baseID: preferred magnitude not found in database: %s",
			               prefix.c_str(), evt->preferredMagnitudeID().c_str());
			return false;
		}

		try {
			magnitude = mag->magnitude().value();
		}
		catch ( ... ) {
			mag = NULL;
		}

		db->loadArrivals(origin.get());

		for ( size_t i = 0; i < origin->arrivalCount(); ++i ) {
			DataModel::Arrival *ar = origin->arrival(i);
			if ( ar->pickID().empty() ) continue;
			DataModel::PickPtr pick = DataModel::Pick::Cast(db->getObject(DataModel::Pick::TypeInfo(), ar->pickID()));
			if ( !pick ) {
				SEISCOMP_WARNING("%s: pick not found in database, ignoring: %s",
				                 prefix.c_str(), ar->pickID().c_str());
				continue;
			}

			SEISCOMP_DEBUG("%s: + pick %s", prefix.c_str(), pick->publicID().c_str());
			EventArrival ea;
			ea.pick = pick;
			try { ea.weight = ar->weight(); }
			catch ( ... ) {}
			ea.phase = ar->phase();

			picks.push_back(ea);
		}
	}
	catch ( ... ) {}

	// In case no event has been read from database we check for a given XML
	// file
	if ( !evt ) {
		try {
			std::string xml = cfg.getString(prefix + "xml");
			IO::XMLArchive ar;
			if ( !ar.open(xml.c_str()) ) {
				SEISCOMP_WARNING("%s.xml: could not open file: %s",
				                 prefix.c_str(), xml.c_str());
				return false;
			}

			DataModel::EventParametersPtr ep;
			bool wasEnabled = DataModel::PublicObject::IsRegistrationEnabled();
			DataModel::PublicObject::SetRegistrationEnabled(false);

			ar >> ep;
			if ( !ep ) {
				SEISCOMP_WARNING("%s.xml: no event parameters: %s",
				                 prefix.c_str(), xml.c_str());
				return false;
			}

			if ( ep->eventCount() == 0 ) {
				SEISCOMP_WARNING("%s.xml: no event found: %s",
				                 prefix.c_str(), xml.c_str());
				return false;
			}

			if ( ep->eventCount() > 1 ) {
				SEISCOMP_WARNING("%s.xml: too many events: %s",
				                 prefix.c_str(), xml.c_str());
				return false;
			}

			evt = ep->event(0);

			origin = findOrigin(ep.get(), evt->preferredOriginID());
			if ( !origin ) {
				SEISCOMP_ERROR("%s.xml: preferred origin not found: %s",
				               prefix.c_str(), evt->preferredOriginID().c_str());
				return false;
			}

			mag = findMagnitude(origin.get(), evt->preferredMagnitudeID());
			if ( !mag ) {
				SEISCOMP_ERROR("%s.xml: preferred magnitude not found: %s",
				               prefix.c_str(), evt->preferredMagnitudeID().c_str());
				return false;
			}

			try {
				magnitude = mag->magnitude().value();
			}
			catch ( ... ) {
				mag = NULL;
			}

			for ( size_t i = 0; i < origin->arrivalCount(); ++i ) {
				DataModel::Arrival *ar = origin->arrival(i);
				if ( ar->pickID().empty() ) continue;
				DataModel::PickPtr pick = findPick(ep.get(), ar->pickID());
				if ( !pick ) {
					SEISCOMP_WARNING("%s: pick not found, ignoring: %s ",
					                 prefix.c_str(), ar->pickID().c_str());
					continue;
				}

				SEISCOMP_DEBUG("%s: + pick %s", prefix.c_str(), pick->publicID().c_str());
				EventArrival ea;
				ea.pick = pick;
				try { ea.weight = ar->weight(); }
				catch ( ... ) {}
				ea.phase = ar->phase();

				picks.push_back(ea);
			}

			DataModel::PublicObject::SetRegistrationEnabled(wasEnabled);
		}
		catch ( ... ) {}
	}

	try {
		if ( !time.fromString(cfg.getString(prefix + "time").c_str(), "%F %T.%f") ) {
			SEISCOMP_WARNING("%s.time: invalid time definition: %s",
			                 prefix.c_str(), cfg.getString(prefix + "time").c_str());
			return false;
		}
	}
	catch ( ... ) {
		if ( !origin ) {
			SEISCOMP_WARNING("%s.time: missing definition", prefix.c_str());
			return false;
		}
	}

	if ( origin ) {
		if ( time.valid() ) {
			if ( origin->time() != time ) {
				SEISCOMP_WARNING("%s: mismatching time and event origin time: %s != %s",
				                 prefix.c_str(), origin->time().value().iso().c_str(),
				                 time.iso().c_str());
				return false;
			}
		}
		else
			time = origin->time().value();
	}

	try { place = cfg.getString(prefix + "place"); }
	catch ( ... ) {}

	try { group = cfg.getString(prefix + "group"); }
	catch ( ... ) {}

	/*
	try {
		offset = cfg.getDouble(prefix + "offset");
	}
	catch ( ... ) {
		offset = 0;
	}

	try {
		duration = cfg.getDouble(prefix + "duration");
		if ( duration <= 0 ) {
			SEISCOMP_WARNING("%s.duration: value out of bounds: %f",
			                 prefix.c_str(), duration);
			return false;
		}
	}
	catch ( ... ) {
		SEISCOMP_WARNING("%s.duration: missing definition", prefix.c_str());
		return false;
	}
	*/

	try {
		magnitude = cfg.getDouble(prefix + "magnitude");
	}
	catch ( ... ) {
		if ( !mag ) {
			SEISCOMP_WARNING("%s.magnitude: missing definition", prefix.c_str());
			return false;
		}
	}

	try {
		deltaM = cfg.getDouble(prefix + "deltaM");
	}
	catch ( ... ) {}

	try {
		latitude = cfg.getDouble(prefix + "latitude");
		if ( latitude < -90 || latitude > 90 ) {
			SEISCOMP_WARNING("%s.latitude: value out of bounds [-90;90]: %f",
			                 prefix.c_str(), latitude);
			return false;
		}
	}
	catch ( ... ) {
		SEISCOMP_WARNING("%s.latitude: missing definition", prefix.c_str());
		return false;
	}

	try {
		longitude = cfg.getDouble(prefix + "longitude");
		if ( longitude < -180 || longitude > 180 ) {
			SEISCOMP_WARNING("%s.longitude: value out of bounds [-180;180]: %f",
			                 prefix.c_str(), longitude);
			return false;
		}
	}
	catch ( ... ) {
		SEISCOMP_WARNING("%s.longitude: missing definition", prefix.c_str());
		return false;
	}

	try {
		depth = cfg.getDouble(prefix + "depth");
		if ( depth < -100 || depth > 800 ) {
			SEISCOMP_WARNING("%s.depth: value out of bounds [-100;800]: %f",
			                 prefix.c_str(), depth);
			return false;
		}
	}
	catch ( ... ) {
		SEISCOMP_WARNING("%s.depth: missing definition", prefix.c_str());
		return false;
	}

	try {
		Environment *env = Environment::Instance();
		dataURI = env->absolutePath(cfg.getString(prefix + "data"));
	}
	catch ( ... ) {
		SEISCOMP_WARNING("%s.data: missing definition", prefix.c_str());
		return false;
	}

	try {
		waveParams.signalBegin = cfg.getDouble(prefix + "signalBegin");
	}
	catch ( ... ) {
		SEISCOMP_WARNING("%s.signalBegin: missing definition", prefix.c_str());
		return false;
	}

	try {
		waveParams.signalEnd = cfg.getDouble(prefix + "signalEnd");
	}
	catch ( ... ) {
		SEISCOMP_WARNING("%s.signalEnd: missing definition", prefix.c_str());
		return false;
	}

	if ( waveParams.signalBegin > waveParams.signalEnd ) {
		SEISCOMP_WARNING("%s.signalWindow: begin (%fs) is after end (%fs)",
		                 prefix.c_str(), waveParams.signalBegin, waveParams.signalEnd);
		return false;
	}

	if ( waveParams.signalBegin == waveParams.signalEnd ) {
		SEISCOMP_WARNING("%s.signalWindow: empty, nothing to do (%fs == %fs)",
		                 prefix.c_str(), waveParams.signalBegin, waveParams.signalEnd);
		return false;
	}

	// Initialize noise windows with default values
	waveParams.noiseBegin = waveParams.signalBegin;
	waveParams.noiseEnd = waveParams.noiseBegin;
	waveParams.noise2Begin = waveParams.signalBegin;
	waveParams.noise2End = waveParams.noise2Begin;

	cfg.getDouble(waveParams.noiseBegin, prefix + "noiseBegin");
	cfg.getDouble(waveParams.noiseEnd, prefix + "noiseEnd");

	if ( waveParams.noiseBegin > waveParams.noiseEnd ) {
		SEISCOMP_WARNING("%s.noiseWindow: begin (%fs) is after end (%fs)",
		                 prefix.c_str(), waveParams.noiseBegin, waveParams.noiseEnd);
		return false;
	}

	if ( !cfg.getDouble(waveParams.noise2Begin, prefix + "noise2Begin") )
		waveParams.noise2Begin = waveParams.noiseBegin;

	if ( !cfg.getDouble(waveParams.noise2End, prefix + "noise2End") )
		waveParams.noise2End = waveParams.noiseEnd;

	if ( waveParams.noise2Begin > waveParams.noise2End ) {
		SEISCOMP_WARNING("%s.noise2Window: begin (%fs) is after end (%fs)",
		                 prefix.c_str(), waveParams.noise2Begin, waveParams.noise2End);
		return false;
	}

	processData = true;
	cfg.getBool(processData, prefix + "processing.enable");

	// Override event specific filter parameters
	cfg.getInt   (waveParams.filterOrder, prefix + "filter.order");
	cfg.getDouble(waveParams.filterFreqs.first, prefix + "filter.loFreq");
	cfg.getDouble(waveParams.filterFreqs.second, prefix + "filter.hiFreq");
	cfg.getBool  (waveParams.filterBandStop, prefix + "filter.bandStop");
	cfg.getBool  (waveParams.envelopeEnable, prefix + "envelope.enable");
	cfg.getInt   (waveParams.envelopeSamplingFrequency, prefix + "envelope.samplingFrequency");
	cfg.getBool  (waveParams.envelopeDownsampleAverage, prefix + "envelope.resampleAverage");
	cfg.getDouble(waveParams.envelopeHiFreq, prefix + "envelope.hiFreq");
	cfg.getBool  (waveParams.envelopeAcausal, prefix + "envelope.acausal");
	cfg.getBool  (waveParams.processingAcausal, prefix + "processing.acausal");
	cfg.getBool  (waveParams.processingLog, prefix + "processing.logarithm");

	if ( waveParams.envelopeEnable ) {
		waveParams.offset = std::min(waveParams.noiseBegin, std::min(waveParams.noise2Begin, waveParams.signalBegin));
		// Need some more data at the beginning
		if ( waveParams.envelopeHiFreq > 0 )
			waveParams.offset -= 1.5 / waveParams.envelopeHiFreq;

		waveParams.duration = std::max(waveParams.noiseEnd, std::max(waveParams.noise2End, waveParams.signalEnd)) - waveParams.offset;
	}
	else {
		waveParams.offset = waveParams.signalBegin;
		waveParams.duration = waveParams.signalEnd - waveParams.offset;
	}

	return true;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::MasterEvent::readData(const Config &cfg,
                                     const StreamList &streams) {
	IO::RecordStreamPtr rs = IO::RecordStream::Open(dataURI.c_str());
	if ( rs == NULL ) {
		SEISCOMP_WARNING("event.%s.data: could not create/open source",
		                 name.c_str());
		return false;
	}

	Core::TimeSpan margin(cfg.detectorMaxTimeLag);
	// Set data acquisition time window from origin time + offset to origin time + offset + duration
	Core::TimeWindow tw(time-margin + Core::TimeSpan(waveParams.offset), time+margin + Core::TimeSpan(waveParams.offset+waveParams.duration));
	rs->setTimeWindow(tw);

	std::map<std::string, double> gains;

	// Add streams
	StreamList::const_iterator it;
	for ( it = streams.begin(); it != streams.end(); ++it ) {
		DataModel::Stream *cha = *it;
		DataModel::SensorLocation *loc = cha->sensorLocation();
		DataModel::Station *sta = loc->station();
		DataModel::Network *net = sta->network();

		gains[net->code() + "." + sta->code() + "." + loc->code() + "." + cha->code()] = cha->gain();

		rs->addStream(net->code(), sta->code(), loc->code(), cha->code());

		std::string staid = net->code() + "." + sta->code();
		if ( stations[staid] == NULL ) stations[staid] = new StationData(loc);
	}

	std::set<std::string> ids;

	// Read records
	IO::RecordInput inp(rs.get(), Array::DOUBLE, Record::DATA_ONLY);
	RecordPtr rec;
	while ( (rec = inp.next()) != NULL ) {
		std::map<std::string, double>::iterator it;
		std::string id = rec->streamID();
		it = gains.find(id);
		// Not in request list. Can happen if a file source is used which does
		// not filter on requested channels
		if ( it == gains.end() ) {
			if ( ids.find(id) == ids.end() ) {
				SEISCOMP_WARNING("[%s] no gain for channel %s: ignoring",
				                 name.c_str(), id.c_str());
				ids.insert(id);
			}

			continue;
		}

		// Store record in event
		Data &d = data[id];
		if ( d.seq == NULL ) {
			SEISCOMP_DEBUG("[%s] adding channel %s", name.c_str(), id.c_str());
			d.seq = new TimeWindowBuffer(tw);
		}
		d.seq->feed(rec.get());
	}

	bool success = true;
	DataSet::iterator sit;

	ids.clear();

	for ( sit = data.begin(); sit != data.end(); ++sit ) {
		if ( sit->second.seq->empty() ) {
			SEISCOMP_ERROR("%s: data available but not inside the requested "
			               "time window of the event", sit->first.c_str());
			success = false;
			sit->second.dispose();
			continue;
		}

		if ( !merge(sit->second.trace, sit->second.seq) ) {
			SEISCOMP_WARNING("%s: data records could not be merged into a single trace",
			                 sit->first.c_str());
			success = false;
		}

		if ( !trim(sit->second.trace, tw) ) {
			SEISCOMP_WARNING("%s: incomplete trace, not enough data for time window %s~%s",
			                 sit->first.c_str(), tw.startTime().iso().c_str(),
			                 tw.endTime().iso().c_str());
			success = false;
		}

		correct(sit->second.trace, gains[sit->second.trace.streamID()]);

		// Delete record sequence
		sit->second.dispose();

		std::string staid = sit->second.trace.networkCode() + "." +
		                    sit->second.trace.stationCode();
		sit->second.stationData = stations[staid];
		ids.insert(staid);
	}

	StationSet::iterator ssit;

	// Trim station set
	for ( ssit = stations.begin(); ssit != stations.end(); ) {
		if ( ids.find(ssit->first) == ids.end() )
			stations.erase(ssit++);
		else
			++ssit;
	}

	return success;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::prepareTrace(Trace &trace, const WaveParams &params, double *pgv) {
	if ( params.envelopeEnable && params.envelopeSamplingFrequency > trace.samplingFrequency() ) {
		SEISCOMP_WARNING("%s: cannot downsample %.1fHz to %.1fHz with respect to the envelope target "
		                 "sampling frequency",
		                 trace.streamID().c_str(),
		                 trace.samplingFrequency(), (double)params.envelopeSamplingFrequency);
		return false;
	}

	DoubleArray *data = DoubleArray::Cast(trace.data());
	double mean = data->mean();

	// Demean data
	int cnt = data->size();
	for ( int i = 0; i < cnt; ++i ) (*data)[i] -= mean;

	if ( params.processingAcausal ) {
		DoubleArrayPtr tmp;
		double Tzpad = 1.5*params.filterOrder;

		if ( params.filterFreqs.first > 0 && params.filterOrder > 0 )
			Tzpad /= params.filterFreqs.first;
		else
			Tzpad = 0.0;

		Tzpad *= 0.5;

		int numZeros = Tzpad * trace.samplingFrequency();

		if ( numZeros > 0 ) {
			DoubleArray zeros(numZeros);
			zeros.fill(0);

			tmp = new DoubleArray(*data);
			data = tmp.get();

			data->prepend(&zeros);
			data->append(&zeros);

			// Taper data
			double taperLength = 1.0;

			if ( taperLength > 0 ) {
				int numTaperSamples = taperLength * trace.samplingFrequency();

				costaper(data->size(), data->typedData(),
				         numZeros, numZeros + numTaperSamples,
				         data->size() - numZeros - numTaperSamples, data->size() - numZeros);
			}
		}

		std::vector<Math::Complex> spec;
		Math::fft(spec, data->size(), data->typedData());
		double df = 1.0/trace.samplingFrequency();

		if ( params.filterOrder > 0 ) {
			if ( params.filterFreqs.first > 0 && params.filterFreqs.second > 0 )
				bandpass(spec, df, params.filterOrder, params.filterFreqs.first, params.filterFreqs.second);
			else if ( params.filterFreqs.first > 0 )
				hipass(spec, df, params.filterOrder, params.filterFreqs.first);
			else if ( params.filterFreqs.second > 0 )
				lopass(spec, df, params.filterOrder, params.filterFreqs.second);
		}

		if ( params.filterBandStop ) {
			// Bandstop at 50Hz
			if ( trace.samplingFrequency() > 102 )
				bandstop(spec, df, 4, 49, 51);

			if ( trace.samplingFrequency() > 202 )
				bandstop(spec, df, 4, 99, 101);
		}

		// Back to time domain
		{
			std::vector<Seiscomp::Math::Complex> tmp_spec(spec);
			Seiscomp::Math::ifft(data->size(),data->typedData(),tmp_spec);
		}

		if ( pgv != NULL ) {
			int idx = find_absmax(data->size(), data->typedData(), 0, data->size(), 0.0);
			*pgv = fabs(data->get(idx));
		}

		if ( params.envelopeEnable ) {
			DoubleArrayPtr in = new DoubleArray(*data);

			if ( params.envelopeHiFreq > 0 ) {
				gausslopass(spec, df, params.envelopeHiFreq);
				std::vector<Seiscomp::Math::Complex> tmp_spec(spec);
				Seiscomp::Math::ifft(in->size(),in->typedData(),tmp_spec);
			}

			hilbert(spec);
			Seiscomp::Math::ifft(data->size(),data->typedData(),spec);
			::magnitude(in->typedData(),data->size(),data->typedData());

			// Copy back data
			if ( numZeros > 0 ) {
				data = DoubleArray::Cast(trace.data());
				memcpy(data->typedData(), tmp->typedData() + numZeros, data->size()*data->elementSize());
			}

			if ( params.processingLog ) log(*data);

			resample(trace, params.envelopeSamplingFrequency, params.envelopeDownsampleAverage);
		}
		else {
			// Copy back data
			if ( numZeros > 0 ) {
				data = DoubleArray::Cast(trace.data());
				memcpy(data->typedData(), tmp->typedData() + numZeros, data->size()*data->elementSize());
			}

			if ( params.processingLog )	log(*data);
		}
	}
	else {
		if ( params.filterOrder > 0 ) {
			if ( params.filterFreqs.first > 0 && params.filterFreqs.second > 0 ) {
				Math::Filtering::IIR::ButterworthHighLowpass<double> bp(params.filterOrder, params.filterFreqs.first, params.filterFreqs.second);
				bp.setSamplingFrequency(trace.samplingFrequency());
				bp.apply(data->size(), data->typedData());
			}
			else if ( params.filterFreqs.first > 0 ) {
				Math::Filtering::IIR::ButterworthHighpass<double> hp(params.filterOrder, params.filterFreqs.first);
				hp.setSamplingFrequency(trace.samplingFrequency());
				hp.apply(data->size(), data->typedData());
			}
			else if ( params.filterFreqs.second > 0 ) {
				Math::Filtering::IIR::ButterworthLowpass<double> lp(params.filterOrder, params.filterFreqs.second);
				lp.setSamplingFrequency(trace.samplingFrequency());
				lp.apply(data->size(), data->typedData());
			}
		}

		if ( params.filterBandStop ) {
			// Bandstop at 50Hz
			if ( trace.samplingFrequency() > 102 ) {
				ButterworthBandstop<double> bs(4,49,51);
				bs.setSamplingFrequency(trace.samplingFrequency());
				bs.apply(data->size(), data->typedData());
			}

			// Bandstop at 100Hz
			if ( trace.samplingFrequency() > 202 ) {
				ButterworthBandstop<double> bs(4,99,101);
				bs.setSamplingFrequency(trace.samplingFrequency());
				bs.apply(data->size(), data->typedData());
			}
		}

		if ( pgv != NULL ) {
			int idx = find_absmax(data->size(), data->typedData(), 0, data->size(), 0.0);
			*pgv = fabs(data->get(idx));
		}

		if ( params.envelopeEnable ) {
			if ( params.envelopeAcausal ) {
				std::vector<Seiscomp::Math::Complex> spec,spec_in;
				Math::fft(spec, data->size(), data->typedData());
				DoubleArrayPtr in = new DoubleArray(*data);

				if ( params.envelopeHiFreq > 0 ) {
					gausslopass(spec, 1.0/trace.samplingFrequency(), params.envelopeHiFreq);
					spec_in = spec;
					Seiscomp::Math::ifft(in->size(),in->typedData(),spec_in);
				}

				hilbert(spec);
				Seiscomp::Math::ifft(data->size(),data->typedData(),spec);
				::magnitude(in->typedData(),data->size(),data->typedData());
			}
			else {
				/*
				if ( cfg.envelopeHiFreq > 0 ) {
					Math::Filtering::IIR::ButterworthHighLowpass<double> lp(cfg.filterOrder, cfg.envelopeHiFreq);
					lp.setSamplingFrequency(trace.samplingFrequency());
					lp.apply(data->size(), data->typedData());
				}
				*/
				Envelope<double> env(1.0 / (params.envelopeHiFreq > 0?params.envelopeHiFreq:trace.samplingFrequency()));
				env.setSamplingFrequency(trace.samplingFrequency());
				env.apply(data->size(), data->typedData());
			}

			if ( params.processingLog ) log(*data);

			resample(trace, params.envelopeSamplingFrequency, params.envelopeDownsampleAverage);
		}
		else if ( params.processingLog )
			log(*data);
	}

	return prepareSignalWindow(trace, params);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::prepareSignalWindow(Trace &trace, const WaveParams &params) {
	DoubleArray *data = DoubleArray::Cast(trace.data());
	int cnt = data->size();
	double mean;

	if ( params.envelopeEnable ) {
		int noiseBegin = (int)((params.noiseBegin-params.offset) * trace.samplingFrequency() + 0.5);
		int noiseEnd = (int)((params.noiseEnd-params.offset) * trace.samplingFrequency() + 0.5);

		int noise2Begin = (int)((params.noise2Begin-params.offset) * trace.samplingFrequency() + 0.5);
		int noise2End = (int)((params.noise2End-params.offset) * trace.samplingFrequency() + 0.5);

		noiseBegin = std::max(noiseBegin, 0);
		noiseEnd = std::min(noiseEnd, cnt);

		noise2Begin = std::max(noise2Begin, 0);
		noise2End = std::min(noise2End, cnt);

		int noiseSamples = noiseEnd - noiseBegin;
		int noise2Samples = noise2End - noise2Begin;

		OPT(double) means[2];

		if ( noiseSamples > 0 ) {
			// Compute mean of noise window
			means[0] = 0.0;
			for ( int i = noiseBegin; i < noiseEnd; ++i ) *means[0] += (*data)[i];
			*means[0] /= noiseSamples;
		}

		if ( noise2Samples > 0 && noise2Begin != noiseBegin && noise2End != noiseEnd ) {
			// Compute mean of noise window 2
			means[1] = 0.0;
			for ( int i = noise2Begin; i < noise2End; ++i ) *means[1] += (*data)[i];
			*means[1] /= noise2Samples;
		}

		if ( means[0] && !means[1] ) {
			mean = *means[0];
			//SEISCOMP_DEBUG("%s: mean = %f", trace.streamID().c_str(), mean);
		}
		else if ( !means[0] && means[1] ) {
			mean = *means[1];
			//SEISCOMP_DEBUG("%s: mean = %f", trace.streamID().c_str(), mean);
		}
		else if ( means[0] && means[1] ) {
			mean = std::min(*means[0], *means[1]);
			//SEISCOMP_DEBUG("%s: min(mean1,mean2) = %f (mean1 = %f, mean2 = %f)",
			//               trace.streamID().c_str(), mean, *means[0], *means[1]);
		}
		else {
			mean = 0;
			//SEISCOMP_DEBUG("%s: mean removal disabled", trace.streamID().c_str());
		}

		// Remove mean
		trace.mean = mean;
		for ( int i = 0; i < cnt; ++i ) (*data)[i] -= mean;
	}

	double noiseOffset = std::min(params.noiseBegin, std::min(params.noise2Begin, params.signalBegin));

	trace.noiseOffset = (int)((noiseOffset-params.offset) * trace.samplingFrequency() + 0.5);
	trace.noiseCnt = (int)((params.signalEnd-params.offset) * trace.samplingFrequency() + 0.5);

	trace.signalOffset = (int)((params.signalBegin-params.offset) * trace.samplingFrequency() + 0.5);
	trace.signalCnt = (int)((params.signalEnd-params.offset) * trace.samplingFrequency() + 0.5);

	if ( trace.signalOffset < 0 ) trace.signalOffset = 0;
	if ( trace.signalCnt > cnt ) trace.signalCnt = cnt;

	if ( trace.noiseOffset < 0 ) trace.noiseOffset = 0;
	if ( trace.noiseCnt > cnt ) trace.noiseCnt = cnt;

	trace.signalCnt -= trace.signalOffset;
	trace.noiseCnt -= trace.noiseOffset;

	return true;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::MasterEvent::prepareData(const Config &cfg) {
	DataSet::iterator it;
	bool success = true;

	for ( it = data.begin(); it != data.end(); ++it ) {
		Trace &trace = it->second.trace;

		if ( processData && !Detector::prepareTrace(trace, waveParams, &it->second.state.pgv) )
			success = false;
		else
			Detector::prepareSignalWindow(trace, waveParams);

		it->second.state.norm = norm(trace.signalCount(), trace.signal(), it->second.state.normSum, 0);
	}

	return success;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::MasterEvent::dumpConfig(std::ostream &os) const {
	os << "[" << name << "]" << std::endl;
	os << " place     : " << place << std::endl;
	os << " group     : " << group << std::endl;
	os << " negative  : " << doNotSendOrigins << std::endl;
	os << " time      : " << time.iso() << std::endl;
	os << " offset    : " << waveParams.offset << "s" << std::endl;
	os << " duration  : " << waveParams.duration << "s" << std::endl;
	os << " magnitude : " << magnitude << std::endl;
	os << " latitude  : " << latitude << std::endl;
	os << " longitude : " << longitude << std::endl;
	os << " depth     : " << depth << "km" << std::endl;
	os << " prepare   : " << processData << std::endl;
	os << " params {" << std::endl;
	os << "  noise    : " << waveParams.noiseBegin << "s ; " << waveParams.noiseEnd << "s" << std::endl;
	os << "  noise2   : " << waveParams.noise2Begin << "s ; " << waveParams.noise2End << "s" << std::endl;
	os << "  signal   : " << waveParams.signalBegin << "s ; " << waveParams.signalEnd << "s" << std::endl;
	os << "  log      : " << waveParams.processingLog << std::endl;
	os << "  acausal  : " << waveParams.processingAcausal << std::endl;
	os << "  filter   : " << waveParams.filterOrder << ","
	                      << waveParams.filterFreqs.first << ","
	                      << waveParams.filterFreqs.second << std::endl;
	os << "  bandstop : " << waveParams.filterBandStop << std::endl;
	os << "  envelope { " << std::endl;
	os << "   enable  : " << waveParams.envelopeEnable << std::endl;
	os << "   SR      : " << waveParams.envelopeSamplingFrequency << std::endl;
	os << "   SR_avg  : " << waveParams.envelopeDownsampleAverage << std::endl;
	os << "   hiFreq  : " << waveParams.envelopeHiFreq << std::endl;
	os << "   acausal : " << waveParams.envelopeAcausal << std::endl;
	os << "  }" << std::endl;
	os << " }" << std::endl;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Detector::Detector(int argc, char **argv) : Application(argc, argv) {
	setLoadStationsEnabled(true);
	setLoadConfigModuleEnabled(false);

	setMessagingEnabled(true);
	setMessagingUsername("magsdetec");
	setPrimaryMessagingGroup("LOCATION");

	NEW_OPT(_config.waveParams.filterOrder, "filter.order");
	NEW_OPT(_config.waveParams.filterFreqs.first, "filter.loFreq");
	NEW_OPT(_config.waveParams.filterFreqs.second, "filter.hiFreq");
	NEW_OPT(_config.waveParams.filterBandStop, "filter.bandStop");

	NEW_OPT(_config.waveParams.envelopeEnable, "envelope.enable");
	NEW_OPT(_config.waveParams.envelopeSamplingFrequency, "envelope.samplingFrequency");
	NEW_OPT(_config.waveParams.envelopeDownsampleAverage, "envelope.resampleAverage");
	NEW_OPT(_config.waveParams.envelopeHiFreq, "envelope.hiFreq");
	NEW_OPT(_config.waveParams.envelopeAcausal, "envelope.acausal");

	NEW_OPT(_config.waveParams.processingAcausal, "processing.acausal");
	NEW_OPT(_config.waveParams.processingLog, "processing.logarithm");

	NEW_OPT(_config.bufferSize, "processing.bufferSize");
	NEW_OPT(_config.processingInterval, "processing.interval");

	NEW_OPT(_config.maximumLatency, "processing.maximumLatency");

	NEW_OPT(_config.normalizationString, "processing.normalization");

	NEW_OPT(_config.detectorThreshold, "detector.threshold");
	NEW_OPT(_config.detectorChannelThreshold, "detector.channelThreshold");
	NEW_OPT(_config.detectorTimeWindow, "detector.window");
	//NEW_OPT(_config.detectorMaxTimeLag, "detector.maxTimeLag");
	NEW_OPT(_config.detectorMinStationRatio, "detector.minimumStationRatio");
	NEW_OPT(_config.detectorMinChannelRatio, "detector.minimumChannelRatio");
	NEW_OPT(_config.detectorProcessingTimeWindow, "detector.minimumProcessingWindow");
	NEW_OPT(_config.detectorPublicationTimeout, "detector.publicationTimeout");
	NEW_OPT(_config.detectorMaximumStepFrequency, "processing.maximumStepFrequency");

	NEW_OPT(_config.channels, "channels");
	NEW_OPT(_config.events, "events");

	NEW_OPT(_config.dumpWaveforms, "output.waveforms.enable");
	NEW_OPT(_config.dumpTriggerWaveforms, "output.waveforms.mseed");
	NEW_OPT(_config.dumpTriggerFit, "output.fit.enable");
	NEW_OPT(_config.dumpWaveformsPath, "output.waveforms.path");
	NEW_OPT(_config.outputEvents, "output.events.file");

	NEW_OPT_CLI(_config.offline, "Mode", "offline", "Do not connect to messaging system", false, true);
	NEW_OPT_CLI(_config.test, "Messaging", "test", "Do not sent messages", false, true);
	NEW_OPT_CLI(_config.dumpFitFunction, "Test", "dump-fit", "Dump continuous fit function to [event.name].fit and [event.name-cha].fit", &_config.dumpFitFunction);
	NEW_OPT_CLI(_config.dumpFitChannels, "Test", "dump-channels", "Dump alls processed channels for the first n time steps to [event.name-cha-time-step].plot. -1 dumps all time steps", &_config.dumpFitChannels);
	NEW_OPT_CLI(_config.showProgress, "Test", "progress", "Print current processing time to stderr", false, true);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Detector::~Detector() {
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::printVersion() {
	std::cout << name() << " application version: "VERSION << std::endl;
	Application::printVersion();
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::validateParameters() {
	if ( !Application::validateParameters() ) return false;

	if ( _config.bufferSize <= 0 ) {
		std::cerr << "Invalid buffer size: value must be positive and non zero"
		          << std::endl;
		return false;
	}

	if ( _config.normalizationString == "total" )
		_config.normalization = NORM_TOTAL;
	else if ( _config.normalizationString == "trace" )
		_config.normalization = NORM_TRACE;
	else {
		std::cerr << "Invalid option '" << _config.normalizationString
		          << "' in normalization" << std::endl;
		return false;
	}

	if ( _config.events.empty() ) {
		std::cerr << "No master events defined, nothing to do" << std::endl;
		return false;
	}

	bool useEventDB = false;
	for ( size_t i = 0; i < _config.events.size(); ++i ) {
		try {
			std::string baseEventID = configGetString("event." + _config.events[i] + ".baseID");
			if ( !baseEventID.empty() ) {
				useEventDB = true;
				break;
			}
		}
		catch ( ... ) {}
	}

	if ( _config.channels.empty() ) {
		std::cerr << "No channels defined, nothing to do" << std::endl;
		return false;
	}

	for ( size_t i = 0; i < _config.channels.size(); ++i ) {
		std::vector<std::string> toks;
		if ( Core::split(toks, _config.channels[i].c_str(), ".", false) != 4 ) {
			std::cerr << "Channel definition at position " << (i+1)
			          << " incomplete, 4 tokens required" << std::endl;
			return false;
		}

		_channels.push_back(DataModel::WaveformStreamID(toks[0], toks[1], toks[2], toks[3], ""));
	}

	Environment *env = Environment::Instance();

	// Convert paths
	_config.outputEvents = env->absolutePath(_config.outputEvents);

	_config.dumpWaveformsPath = env->absolutePath(_config.dumpWaveformsPath);
	if ( !_config.dumpWaveformsPath.empty() && *_config.dumpWaveformsPath.rbegin() != '/' )
		_config.dumpWaveformsPath += '/';

	if ( _config.dumpWaveforms ) {
		if ( !Util::pathExists(_config.dumpWaveformsPath) ) {
			if ( !Util::createPath(_config.dumpWaveformsPath) ) {
				SEISCOMP_ERROR("Unable to create waveform output directory: %s",
				               _config.dumpWaveformsPath.c_str());
				return false;
			}
		}
	}

	if ( !isInventoryDatabaseEnabled() && !useEventDB )
		setDatabaseEnabled(false, false);

	if ( _config.offline ) {
		_config.test = true;
		setMessagingEnabled(false);
	}

	return true;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::init() {
	if ( !Application::init() ) return false;

	_logOrigins = addOutputObjectLog("origin", primaryMessagingGroup());

	for ( size_t i = 0; i < _config.events.size(); ++i ) {
		MasterEventPtr evt = new MasterEvent;
		evt->name = _config.events[i];

		// Copy default processing parameters
		evt->waveParams = _config.waveParams;

		if ( !evt->read(configuration(), query(), "event." + evt->name + ".") ) {
			SEISCOMP_ERROR("Incomplete master event definition for '%s'",
			               evt->name.c_str());
			_masterEvents.clear();
			return false;
		}

		evt->dumpConfig(std::cerr);

		StreamList streams;
		resolveStreams(streams, evt->time);
		if ( streams.empty() ) {
			SEISCOMP_ERROR("No channels found in inventory to read master event '%s'",
			               evt->name.c_str());
			_masterEvents.clear();
			return false;
		}

		SEISCOMP_INFO("Reading master event %s", evt->name.c_str());
		if ( !evt->readData(_config, streams) ) {
			SEISCOMP_ERROR("Failed to read data for master event '%s'",
			               evt->name.c_str());
			_masterEvents.clear();
			return false;
		}

		if ( evt->data.empty() ) {
			SEISCOMP_WARNING("No data read for master event '%s': ignoring",
			                 evt->name.c_str());
			continue;
		}

		if ( !evt->prepareData(_config) ) {
			SEISCOMP_ERROR("Failed to prepare data for master event '%s'",
			               evt->name.c_str());
			_masterEvents.clear();
			return false;
		}

		evt->processing.reset();

		_masterEvents.push_back(evt);
		if ( !evt->group.empty() )
			_eventGroups[evt->group].events.push_back(evt);

		// Find minimum required event duration used later for real time
		// processing
		if ( _masterEvents.size() == 1 )
			_minimumRequiredWindowLength = evt->waveParams.duration;
		else
			_minimumRequiredWindowLength = std::min(_minimumRequiredWindowLength, evt->waveParams.duration);
	}

	StreamList streams;
	resolveStreams(streams, Core::Time::GMT());

	if ( streams.empty() ) {
		SEISCOMP_ERROR("Could not match any requested stream with current inventory, nothing to do");
		_masterEvents.clear();
		return false;
	}

	// Subscribe to streams
	StreamList::const_iterator it;
	for ( it = streams.begin(); it != streams.end(); ++it ) {
		DataModel::Stream *cha = *it;
		DataModel::SensorLocation *loc = cha->sensorLocation();
		DataModel::Station *sta = loc->station();
		DataModel::Network *net = sta->network();
		addStream(net->code(), sta->code(), loc->code(), cha->code());

		RTTracePtr trace = new RTTrace;
		trace->setGain(1.0);

		std::string id = net->code()+"."+sta->code()+"."+loc->code()+"."+cha->code();
		_waveforms[id] = trace;
	}

	MasterEvents::iterator mit;
	for ( mit = _masterEvents.begin(); mit != _masterEvents.end(); ++mit ) {
		const MasterEvent &evt = **mit;
		DataSet::const_iterator dit;

		for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
			dit->second.dataState.source = NULL;

			// Do we have waveforms for this master event trace?
			Waveforms::const_iterator wit = _waveforms.find(dit->first);
			if ( wit == _waveforms.end() ) continue;

			dit->second.dataState.source = wit->second.get();
		}

		if ( _config.dumpFitFunction )
			evt.processing.fitOutput.open((evt.name + ".fit").c_str());
	}

	_recordsReceived = false;

	// Enable timer triggered each second
	enableTimer(1);

	_currentProcessingTicker = _config.processingInterval;

	return true;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::run() {
	return Application::run();
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::done() {
	// Flush processing
	process(0.0, 0.0);
	flushResults(true);

	Application::done();

	//dumpWaveformState();
	if ( _config.dumpWaveforms ) dumpMasterEvents();

	// Write out index files
	if ( _config.dumpFitFunction ) {
		std::ofstream ofs("events.list");
		MasterEvents::iterator mit;
		for ( mit = _masterEvents.begin(); mit != _masterEvents.end(); ++mit ) {
			const MasterEvent &evt = **mit;
			if ( evt.processing.fitOutput.is_open() ) {
				ofs << evt.name << std::endl;

				std::ofstream cha_ofs((std::string("events-") + evt.name + ".list").c_str());
				OutputMap::iterator oit;
				for ( oit = evt.processing.fitChannelOutputs.begin();
				      oit != evt.processing.fitChannelOutputs.end(); ++oit ) {
					cha_ofs << oit->first << std::endl;
				}
			}
		}
		ofs.close();
	}
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::dumpWaveformState() {
	Waveforms::iterator it;
	for ( it = _waveforms.begin(); it != _waveforms.end(); ++it ) {
		std::cerr << it->first << " " << it->second->startTime().toString("%F %T.%6f")
		          << " " << it->second->endTime().toString("%F %T.%6f") << " "
		          << it->second->sampleCount() << std::endl;
	}
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::dumpMasterEvents() {
	MasterEvents::iterator it;

	std::string path = _config.dumpWaveformsPath;

	for ( it = _masterEvents.begin(); it != _masterEvents.end(); ++it ) {
		const MasterEvent &evt = **it;
		SEISCOMP_INFO("Dumping master event %s to %s%s.mseed",
		              evt.name.c_str(), path.c_str(), evt.name.c_str());

		std::ofstream ofs;
		ofs.open((path + evt.name + ".mseed").c_str());

		DataSet::const_iterator sit;
		for ( sit = evt.data.begin(); sit != evt.data.end(); ++sit ) {
			const Trace &trace = sit->second.trace;
			IO::MSeedRecord(trace).write(ofs);
		}

		ofs.close();
	}
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::dumpWaveforms(const MasterEvent &evt, const Core::Time &time,
                             double fit) {
	DataSet::const_iterator it;
	std::ofstream ofs;
	char buf[16];

	std::string path = _config.dumpWaveformsPath + evt.name;
	if ( !Util::pathExists(path) ) {
		if ( !Util::createPath(path) ) {
			SEISCOMP_ERROR("Unable to create waveform output directory: %s",
			               path.c_str());
			return;
		}
	}

	path += '/';

	sprintf(buf, "%.4f", fit*100.0);
	ofs.open((path + time.toString("%Y%m%d%H%M%S%6f") + "-" + buf + ".mseed").c_str());

	for ( it = evt.data.begin(); it != evt.data.end(); ++it ) {
		if ( !it->second.dataState.trace ) continue;
		Trace &trace = *it->second.dataState.trace;
		trace.setNetworkCode(it->second.trace.networkCode());
		trace.setStationCode(it->second.trace.stationCode());
		// Set to "continuous trace"
		trace.setLocationCode("CT");
		trace.setChannelCode(it->second.trace.channelCode());
		IO::MSeedRecord(trace).write(ofs);
	}

	ofs.close();
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::handleRecord(Record *rec) {
	RecordPtr tmp(rec);

	std::string id = rec->streamID();
	Waveforms::iterator rit = _waveforms.find(id);
	// Ignore record if stream is not registered
	if ( rit == _waveforms.end() ) return;

	RTTracePtr trace = rit->second;

	// Initialize trace if necessary
	if ( !trace->startTime().valid() ) {
		trace->setSamplingFrequency(rec->samplingFrequency());
		trace->setTimeSpan(Core::TimeSpan(_config.bufferSize,0));
	}

	trace->feed(rec);

	if ( _config.processingInterval > 0 )
		_recordsReceived = true;
	else
		process(_config.detectorProcessingTimeWindow,
		        _config.maximumLatency);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::handleTimeout() {
	// Check for processing timeout
	if ( _currentProcessingTicker > 0 ) {
		--_currentProcessingTicker;
		if ( _currentProcessingTicker ) {
			_currentProcessingTicker = _config.processingInterval;

			if ( _recordsReceived ) {
				process(_config.detectorProcessingTimeWindow,
				        _config.maximumLatency);
				_recordsReceived = false;
			}
		}
	}

	// Check for queued result timeouts
	if ( _config.detectorPublicationTimeout >= 0 )
		flushResults(false);
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::flushResults(bool force) {
	// Check for pending events
	EventGroups::iterator it;
	for ( it = _eventGroups.begin(); it != _eventGroups.end(); ++it ) {
		EventGroup &group = it->second;

		ResultQueue::iterator rit;
		for ( rit = group.queuedResults.begin(); rit != group.queuedResults.end(); ) {
			ProcessingResult &res = *rit;
			--res.timeout;

			// Publish event
			if ( force || res.timeout <= 0 ) {
				ProcessingResult other;
				ProcessingResult best;
				bool foundFirst = false;

				SEISCOMP_DEBUG("Event %s in group %s timed out -> publish best of this group",
				               res.event->name.c_str(), it->first.c_str());

				while ( group.popPending(other, res.time, _config.detectorTimeWindow) ) {
					if ( foundFirst ) {
						best = other;
						foundFirst = false;
					}
					else if ( other.fit > best.fit )
						best = other;
				}

				publish(&best);
				SEISCOMP_DEBUG("Events still pending in group %s: %d",
				               it->first.c_str(),
				               (int)group.queuedResults.size());
				// Reset iterator
				rit = group.queuedResults.begin();
			}
			else
				++rit;
		}
	}
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::process(double minTimeWindow, double latency) {
	MasterEvents::iterator eit;
	for ( eit = _masterEvents.begin(); eit != _masterEvents.end(); ++eit ) {
		MasterEvent &evt = **eit;

		// Check for unprocessed common blocks
		Core::TimeWindow  tw;
		Core::Time        minStart, maxEnd;
		DataSet::iterator dit;
		bool first = true;
		bool ignoreEmptyTraces = true;
		int  channelCount = 0;

		if ( latency > 0.0 ) {
			// Find maximum end time of record
			for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
				// Ignore empty traces
				RTTrace *trace = dit->second.dataState.source;
				if ( !trace->startTime().valid() ) continue;

				if ( !minStart.valid() ) minStart = trace->startTime();
				else if ( minStart > trace->startTime() ) minStart = trace->startTime();
				if ( !maxEnd.valid() ) maxEnd = trace->endTime();
				else if ( maxEnd < trace->endTime() ) maxEnd = trace->endTime();
			}

			ignoreEmptyTraces = double(maxEnd-minStart) > latency;
		}

		for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
			// Ignore empty traces
			RTTrace *trace = dit->second.dataState.source;
			if ( !trace->startTime().valid() ) {
				// If empty traces are not to be ignored, skip processing for
				// this event
				if ( !ignoreEmptyTraces ) {
					channelCount = 0;
					break;
				}

				continue;
			}

			// Ignore data with too high latency
			if ( latency > 0 && (double)(maxEnd-trace->endTime()) > latency )
				continue;

			if ( first ) {
				tw.setStartTime(trace->startTime());
				tw.setEndTime(trace->endTime());
				first = false;
			}
			else
				tw = tw.overlap(Core::TimeWindow(trace->startTime(), trace->endTime()));

			++channelCount;
		}

		// No channels
		if ( !channelCount ) continue;

		if ( evt.processedTimeWindow.endTime().valid() ) {
			Core::Time contStartTime = evt.processedTimeWindow.endTime() - Core::TimeSpan(evt.waveParams.duration);
			if ( tw.startTime() < contStartTime ) tw.setStartTime(contStartTime);
			// Processing gap, what to do?
			else if ( tw.startTime() > contStartTime ) {
				SEISCOMP_WARNING("Gap in fitting function for master event %s: %s %s: reset processing",
				                 evt.name.c_str(),
				                 contStartTime.toString("%F %T.%6f").c_str(),
				                 tw.startTime().toString("%F %T.%6f").c_str());
				evt.processing.reset();
			}
		}

		// Nothing to do, time window too short
		if ( tw.length() < std::max(minTimeWindow, evt.waveParams.duration) ) continue;

		// Nothing to do
		if ( tw.endTime() <= evt.processedTimeWindow.endTime() ) continue;

		// Merge processing time window with current time window
		evt.processedTimeWindow = evt.processedTimeWindow | tw;

		SEISCOMP_DEBUG("Processing event %s with %d/%d channels and time window %s %s",
		               evt.name.c_str(), channelCount, (int)evt.data.size(),
		               tw.startTime().toString("%F %T.%6f").c_str(),
		               tw.endTime().toString("%F %T.%6f").c_str());
		process(evt, tw);

		SEISCOMP_DEBUG("Processed data for event %s from %s to %s",
		               evt.name.c_str(),
		               evt.processedTimeWindow.startTime().toString("%F %T.%6f").c_str(),
		               evt.processedTimeWindow.endTime().toString("%F %T.%6f").c_str());
	}

	// Check for pending events which can be flushed
	EventGroups::iterator git;
	for ( git = _eventGroups.begin(); git != _eventGroups.end(); ++git ) {
		EventGroup &group = git->second;
		ResultQueue::const_iterator rit;
		// Iterate over all pending events
		for ( rit = group.queuedResults.begin(); rit != group.queuedResults.end(); ) {
			const ProcessingResult &res = *rit;

			int outstanding = 0;
			MasterEvents::const_iterator mit;
			const ProcessingResult *bestGroupEvent = &res;
			double bestGroupFit = res.fit;

			for ( mit = group.events.begin(); mit != group.events.end(); ++mit ) {
				MasterEvent *master = mit->get();
				if ( master == res.event ) continue;

				const ProcessingResult *res2 = group.pending(master, res.time, _config.detectorTimeWindow);
				if ( res2 != NULL ) {
					// Check for maximum fit
					if ( res2->fit > bestGroupFit ) {
						bestGroupEvent = res2;
						bestGroupFit = res2->fit;
					}
				}
				else if ( master->processedTimeWindow.endTime() < res.time+Core::TimeSpan(_config.detectorTimeWindow+master->waveParams.offset+master->waveParams.duration) ) {
					++outstanding;
					continue;
				}
			}

			if ( !outstanding ) {
				SEISCOMP_DEBUG("No outstanding event processing for group %s, "
				               "publish event with highest fit", git->first.c_str());
				publish(bestGroupEvent);
				group.removePending(res.time, _config.detectorTimeWindow);
				// Reset iterator
				rit = group.queuedResults.begin();
			}
			else
				++rit;
		}
	}
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::resolveStreams(StreamList &streams, const Core::Time &time) {
	DataModel::Inventory *inv = Client::Inventory::Instance()->inventory();

	WaveformIDList::iterator it;
	for ( it = _channels.begin(); it != _channels.end(); ++it ) {
		const DataModel::WaveformStreamID &wid = *it;

		size_t cnt = streams.size();

		for ( size_t n = 0; n < inv->networkCount(); ++n ) {
			DataModel::Network *net = inv->network(n);
			if ( net->start() > time ) continue;
			try { if ( net->end() < time ) continue; }
			catch ( ... ) {}

			if ( !Core::wildcmp(wid.networkCode(), net->code()) ) continue;

			for ( size_t s = 0; s < net->stationCount(); ++s ) {
				DataModel::Station *sta = net->station(s);
				if ( sta->start() > time ) continue;
				try { if ( sta->end() < time ) continue; }
				catch ( ... ) {}

				if ( !Core::wildcmp(wid.stationCode(), sta->code()) ) continue;

				for ( size_t l = 0; l < sta->sensorLocationCount(); ++l ) {
					DataModel::SensorLocation *loc = sta->sensorLocation(l);
					if ( loc->start() > time ) continue;
					try { if ( loc->end() < time ) continue; }
					catch ( ... ) {}

					if ( !Core::wildcmp(wid.locationCode(), loc->code()) ) continue;

					// Use all three components
					if ( wid.channelCode().length() == 2 ) {
						DataModel::ThreeComponents tc;
						DataModel::getThreeComponents(tc, loc, wid.channelCode().c_str(), time);
						if ( tc.comps[DataModel::ThreeComponents::Vertical] == NULL ) {
							SEISCOMP_WARNING("%s.%s.%s.%s: vertical component not found in inventory",
							                 net->code().c_str(), sta->code().c_str(),
							                 loc->code().c_str(), wid.channelCode().c_str());
							continue;
						}
						else if ( tc.comps[DataModel::ThreeComponents::FirstHorizontal] == NULL ) {
							SEISCOMP_WARNING("%s.%s.%s.%s: first horizontal component not found in inventory",
							                 net->code().c_str(), sta->code().c_str(),
							                 loc->code().c_str(), wid.channelCode().c_str());
							continue;
						}
						else if ( tc.comps[DataModel::ThreeComponents::SecondHorizontal] == NULL ) {
							SEISCOMP_WARNING("%s.%s.%s.%s: second horizontal component not found in inventory",
							                 net->code().c_str(), sta->code().c_str(),
							                 loc->code().c_str(), wid.channelCode().c_str());
							continue;
						}
						else {
							for ( int i = 0; i < 3; ++i ) {
								if ( !checkGain(tc.comps[i]) )
									SEISCOMP_WARNING("%s.%s.%s.%s: invalid gain, ignoring channel",
									                 net->code().c_str(), sta->code().c_str(),
									                 loc->code().c_str(), tc.comps[i]->code().c_str());
								else
									streams.push_back(tc.comps[i]);
							}
						}
					}
					else {
						for ( size_t s = 0; s < loc->streamCount(); ++s ) {
							DataModel::Stream *cha = loc->stream(s);
							if ( cha->start() > time ) continue;
							try { if ( cha->end() < time ) continue; }
							catch ( ... ) {}

							if ( !Core::wildcmp(wid.channelCode(), cha->code()) ) continue;

							if ( !checkGain(cha) ) {
								SEISCOMP_WARNING("%s.%s.%s.%s: invalid gain, ignoring channel",
								                 net->code().c_str(), sta->code().c_str(),
								                 loc->code().c_str(), cha->code().c_str());
								continue;
							}

							streams.push_back(cha);
						}
					}
				}
			}
		}

		if ( streams.size() <= cnt ) {
			SEISCOMP_WARNING("No matching streams in inventory for %s.%s.%s.%s",
			                 wid.networkCode().c_str(), wid.stationCode().c_str(),
			                 wid.locationCode().c_str(), wid.channelCode().c_str());
		}
	}
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::processWaveforms() {
	Core::TimeWindow tw;

	Waveforms::iterator it;
	bool first = true;
	for ( it = _waveforms.begin(); it != _waveforms.end(); ++it ) {
		// Ignore empty traces
		if ( !it->second->startTime().valid() ) continue;

		if ( first ) {
			tw.setStartTime(it->second->startTime());
			tw.setEndTime(it->second->endTime());
			first = false;
		}
		else
			tw = tw.overlap(Core::TimeWindow(it->second->startTime(), it->second->endTime()));
	}

	SEISCOMP_INFO("Processing common time window %s %s",
	              tw.startTime().toString("%F %T.%6f").c_str(),
	              tw.endTime().toString("%F %T.%6f").c_str());

	MasterEvents::iterator eit;
	for ( eit = _masterEvents.begin(); eit != _masterEvents.end(); ++eit ) {
		const MasterEvent &evt = **eit;
		SEISCOMP_DEBUG(" Processing master event %s", evt.name.c_str());
		process(evt, tw);
		if ( isExitRequested() ) break;
	}
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::process(const MasterEvent &evt, Core::TimeWindow tw) {
	//tw.setStartTime(Core::Time::FromString("2012-11-12 11:15:00.0", "%F %T.%f"));
	//tw.setStartTime(Core::Time::FromString("2013-02-17 20:07:10.00", "%F %T.%f"));
	//tw.setEndTime(Core::Time::FromString("2013-02-17 20:08:22.00", "%F %T.%f"));

	if ( tw.length() < evt.waveParams.duration ) {
		SEISCOMP_DEBUG("  time window is too short for event duration: %.2f < %.2f",
		               tw.length(), evt.waveParams.duration);
		return;
	}

	double lowestSampleFrequency = -1;

	int stationCount = 0;
	int channelCount = 0;

	DataSet::const_iterator dit;
	StationSet::const_iterator sit;

	for ( sit = evt.stations.begin(); sit != evt.stations.end(); ++sit )
		sit->second->channelsUsed = 0;

	for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
		dit->second.state.isUsed = false;

		if ( dit->second.dataState.source == NULL ) {
			SEISCOMP_DEBUG("  %s: no waveforms: ignoring", dit->first.c_str());
			continue;
		}

		if ( !dit->second.dataState.source->startTime().valid() ) {
			SEISCOMP_DEBUG("  %s: empty waveforms: ignoring", dit->first.c_str());
			continue;
		}

		double gain;

		DataModel::Stream *stream = Client::Inventory::Instance()->getStream(
		                                dit->second.trace.networkCode(),
		                                dit->second.trace.stationCode(),
		                                dit->second.trace.locationCode(),
		                                dit->second.trace.channelCode(),
		                                tw.startTime());

		if ( stream == NULL ) {
			SEISCOMP_WARNING("  %s: stream not found in inventory for epoch %s",
			                 dit->first.c_str(), tw.startTime().toString("%F %T").c_str());
			continue;
		}

		try {
			gain = stream->gain();
		}
		catch ( ... ) {
			SEISCOMP_WARNING("  %s: no gain for epoch %s found",
			                 dit->first.c_str(), tw.startTime().toString("%F %T").c_str());
			continue;
		}

		dit->second.dataState.gain = gain;
		dit->second.dataState.trace = new Trace();
		dit->second.state.isUsed = true;
		++dit->second.stationData->channelsUsed;
		++channelCount;

		if ( lowestSampleFrequency < 0 )
			lowestSampleFrequency = dit->second.dataState.source->samplingFrequency();
		else
			lowestSampleFrequency = std::min(lowestSampleFrequency, dit->second.dataState.source->samplingFrequency());
	}

	if ( evt.waveParams.envelopeEnable && evt.waveParams.envelopeSamplingFrequency > 0 )
		lowestSampleFrequency = evt.waveParams.envelopeSamplingFrequency;

	if ( _config.detectorMaximumStepFrequency > 0 && lowestSampleFrequency > _config.detectorMaximumStepFrequency )
		lowestSampleFrequency = _config.detectorMaximumStepFrequency;

	if ( !channelCount ) {
		SEISCOMP_DEBUG("  no channel data available");
		return;
	}

	for ( sit = evt.stations.begin(); sit != evt.stations.end(); ++sit )
		if ( sit->second->channelsUsed > 0 ) ++stationCount;

	if ( stationCount < ((int)evt.stations.size()*_config.detectorMinStationRatio+50)/100 ) {
		SEISCOMP_WARNING("  not enough stations available to process: %d/%d",
		                 stationCount, (int)evt.stations.size());
		return;
	}

	if ( channelCount < ((int)evt.data.size()*_config.detectorMinChannelRatio+50)/100 ) {
		SEISCOMP_WARNING("  not enough channels available to process: %d/%d",
		                 channelCount, (int)evt.data.size());
		return;
	}

	Core::Time time = tw.startTime();
	Core::Time maxTime = tw.endTime() - Core::TimeSpan(evt.waveParams.duration);
	Core::TimeSpan step;

	if ( _config.detectorMaxTimeLag > 0 ) {
		step = Core::TimeSpan(_config.detectorMaxTimeLag*0.9);
		lowestSampleFrequency = 1.0 / (double)step;
	}
	else
		step = Core::TimeSpan(1.0/lowestSampleFrequency);

	SEISCOMP_DEBUG("  Time step is %d.%06ds", (int)step.seconds(), (int)step.microseconds());

	long   lastsecs = 0;
	int    searchMaximumSteps = (int)(_config.detectorTimeWindow * lowestSampleFrequency);
	Core::TimeSpan margin(_config.detectorMaxTimeLag);

	while ( time < maxTime && !isExitRequested() ) {
		if ( lastsecs != time.seconds()/10 ) {
			if ( _config.showProgress )
				std::cerr << evt.name << "    " << time.toString("%F %T.%1f") << std::endl;
			SEISCOMP_DEBUG("  time is %s", time.toString("%F %T").c_str());
			lastsecs = time.seconds()/10;
		}

		FitNorm globalFit;

		channelCount = 0;

		// Pass 0: Reset used channels items
		for ( sit = evt.stations.begin(); sit != evt.stations.end(); ++sit )
			sit->second->channelsUsed = 0;

		// Pass 1: Cut and prepare data
		for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
			if ( !dit->second.dataState.trace ) continue;

			Core::TimeWindow eventTw(time-margin, time+margin + Core::TimeSpan(evt.waveParams.duration));

			const RTTrace &rt_data = *dit->second.dataState.source;

			// Extract data
			if ( !rt_data.extract(*dit->second.dataState.trace, eventTw) ) {
				SEISCOMP_WARNING("  %s: data extraction failed", dit->first.c_str());
				dit->second.state.isUsed = false;
				continue;
			}

			correct(*dit->second.dataState.trace, dit->second.dataState.gain);
			prepareTrace(*dit->second.dataState.trace, evt.waveParams, &dit->second.dataState.pgv);

			if ( dit->second.trace.samplingFrequency() != dit->second.dataState.trace->samplingFrequency() ) {
				SEISCOMP_WARNING("  %s: incompatible sampling frequencies: master(%.1f) != data(%.1f)",
				                 dit->first.c_str(),
				                 dit->second.trace.samplingFrequency(),
				                 dit->second.dataState.trace->samplingFrequency());
				dit->second.dataState.trace = NULL;
				dit->second.state.isUsed = false;
				continue;
			}

			dit->second.state.isUsed = true;
			++channelCount;
		}

		int fitStationCount = 0;
		int fitChannelCount = 0;

		// Pass 2: compute fit and characteristic function
		for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
			if ( !dit->second.state.isUsed ) continue;

			Trace &trace = *dit->second.dataState.trace;

			int cnt = std::min(trace.signalCount(), dit->second.trace.signalCount());
			const double *data_samples = trace.signal();
			const double *master_samples = dit->second.trace.signal();

			//int maxTimeLag = (int)(_config.detectorMaxTimeLag * dit->second.trace.samplingFrequency());

			FitNorm fitNorm;

			dit->second.state.lag = crosscorr(master_samples, cnt, data_samples, 2, fitNorm);
			//dit->second.state.lag = 0;
			//goodnessOfFit(master_samples, cnt, data_samples, fitNorm);

			dit->second.dataState.norm = norm(cnt, data_samples, dit->second.dataState.normSum, dit->second.state.lag);
			dit->second.state.norm = norm(cnt, master_samples, dit->second.state.normSum, -dit->second.state.lag);

			dit->second.state.fitRaw = fitNorm();

			FitNorm traceNorm = FitNorm(fitNorm).scale(1.0 / (dit->second.dataState.norm*dit->second.state.norm));
			dit->second.state.traceFit = traceNorm();
			if ( dit->second.state.traceFit > 1.0 ) dit->second.state.traceFit = 1.0;

			if ( _config.dumpFitChannels < 0 || evt.processing.loop < (uint64_t)_config.dumpFitChannels ) {
				cnt = std::min(trace.noiseCount(), dit->second.trace.noiseCount());
				data_samples = trace.noise();
				master_samples = dit->second.trace.noise();

				int d_ofs = dit->second.state.lag;
				int m_ofs = 0;
				if ( d_ofs < 0 ) {
					m_ofs = -d_ofs;
					d_ofs = 0;
					cnt -= m_ofs;
				}
				else if ( d_ofs > 0 )
					cnt -= d_ofs;

				double scaleD = 1.0 / dit->second.dataState.norm;
				double scaleM = 1.0 / dit->second.state.norm;

				std::ofstream ofs;
				char buf[128];
				snprintf(buf, 127, "%08ld", evt.processing.loop);
				ofs.open((evt.name + "-" + dit->first + "-" + buf + "-" + time.toString("%Y%m%d%H%M%S%6f") + ".plot").c_str());

				ofs << "# mean = " << dit->second.trace.mean << std::endl;
				ofs << "# master_samples = " << dit->second.trace.signalCount() << std::endl;

				Core::Time t = time;
				double dt = 1.0 / trace.samplingFrequency();

				for ( int i = 0; i < cnt; ++i, t += dt )
					ofs << scaleD*data_samples[i+d_ofs]
					    << " " << scaleM*master_samples[i+m_ofs]
					    << " " << t.toString("%Y%m%d%H%M%S%6f")
					    << std::endl;
				ofs.close();
			}

			dit->second.state.fit = dit->second.state.traceFit;

			if ( dit->second.state.traceFit < _config.detectorChannelThreshold ) {
				dit->second.state.isUsed = false;
				continue;
			}

			//dit->second.state.timeLag = dit->second.state.lag / dit->second.trace.samplingFrequency();
			dit->second.state.timeLag = 0;

			if ( _config.normalization == NORM_TRACE )
				globalFit += traceNorm;
			else
				globalFit += fitNorm;

			++fitChannelCount;
			++dit->second.stationData->channelsUsed;
		}

		// Count used stations
		for ( sit = evt.stations.begin(); sit != evt.stations.end(); ++sit )
			if ( sit->second->channelsUsed > 0 ) ++fitStationCount;

		double fit;
		double totalScale = 1.0;
		double traceScale = 1.0;

		if ( fitStationCount > 0 ) {
			if ( _config.normalization == NORM_TRACE ) {
				traceScale = 1.0 / fitChannelCount;
				globalFit.scale(traceScale);
			}
			else {
				double total_norm = 0.0;
				double m_total_norm = 0.0;

				for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
					if ( !dit->second.state.isUsed ) continue;
					total_norm += dit->second.dataState.normSum;
					m_total_norm += dit->second.state.normSum;
				}

				totalScale = 1.0 / (sqrt(total_norm)*sqrt(m_total_norm));
				globalFit.scale(totalScale);
			}

			fit = globalFit();
			if ( fit > 1.0 ) fit = 1.0;

			/** Check for good station/channel coverage **/

			// We need at least X% of the master stations
			if ( fitStationCount < ((int)evt.stations.size()*_config.detectorMinStationRatio+50)/100 )
				fit = 0;
			// We need at least Y% of available channels fitting
			else if ( fitChannelCount < ((int)evt.data.size()*_config.detectorMinChannelRatio+50)/100 )
				fit = 0;
		}
		else
			fit = 0;

		if ( fit >= _config.detectorThreshold && evt.processing.triggerCounter <= 0 ) {
			evt.processing.triggerCounter = searchMaximumSteps;

			if ( _config.dumpTriggerFit ) {
				std::string path = _config.dumpWaveformsPath + evt.name;
				bool error = false;
				if ( !Util::pathExists(path) ) {
					if ( !Util::createPath(path) ) {
						SEISCOMP_ERROR("Unable to create waveform plot output directory: %s",
						               path.c_str());
						error = true;
					}
				}

				if ( !error ) {
					path += '/';
					evt.processing.triggerFitOutput.open((path + time.toString("%Y%m%d%H%M%S%6f") + "-event.fit").c_str());
				}

				for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
					if ( !dit->second.dataState.trace ) continue;

					OutputMap::iterator oit = evt.processing.triggerFitChannelOutputs.find(dit->first);
					ofstream_ptr ofs;
					if ( oit == evt.processing.triggerFitChannelOutputs.end() ) {
						ofs = ofstream_ptr(new std::ofstream);
						evt.processing.triggerFitChannelOutputs[dit->first] = ofs;
						ofs->open((path + time.toString("%Y%m%d%H%M%S%6f") + "-" + dit->first + ".fit").c_str());
					}
				}
			}
		}

		Core::Time originTime = time + Core::TimeSpan(-evt.waveParams.offset);

		if ( evt.processing.triggerCounter > 0 ) {
			//std::cerr << " " << (fit*100.0) << "%  " << time.toString("%F %T.%f") << "  " << fits << " channels are fitting" << std::endl;
			if ( fit > evt.processing.bestFit ) {
				evt.processing.bestFitChannels.clear();
				evt.processing.bestFittingTime = originTime;
				evt.processing.bestFit = fit;

				int n = 0;
				double sumratio = 0.0;
				for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
					dit->second.dataState.isUsed = dit->second.state.isUsed;
					if ( !dit->second.state.isUsed ) continue;

					double magOfs = log10(dit->second.dataState.pgv/dit->second.state.pgv);
					sumratio += magOfs;
					dit->second.dataState.fit = dit->second.state.fit;
					dit->second.dataState.traceFit = dit->second.state.traceFit;
					dit->second.dataState.timeLag = dit->second.state.timeLag;
					dit->second.dataState.lag = dit->second.state.lag;
					++n;

					evt.processing.bestFitChannels.push_back(ResultItem(&dit->second, dit->second.dataState.fit, evt.magnitude + magOfs));
				}

				evt.processing.bestFitMag = evt.magnitude + sumratio / n;
				if ( evt.deltaM )
					evt.processing.bestFitMag += *evt.deltaM;
				evt.processing.bestFitStationCount = fitStationCount;
				evt.processing.bestFitChannelCount = fitChannelCount;
			}

			if ( _config.dumpTriggerWaveforms ) dumpWaveforms(evt, time, fit);
			if ( _config.dumpTriggerFit ) {
				evt.processing.triggerFitOutput << originTime.iso() << " " << fit << std::endl;

				for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
					if ( !dit->second.dataState.trace ) continue;

					OutputMap::iterator oit = evt.processing.triggerFitChannelOutputs.find(dit->first);
					ofstream_ptr ofs;
					if ( oit == evt.processing.triggerFitChannelOutputs.end() ) continue;
					ofs = oit->second;
					if ( !ofs ) continue;

					if ( dit->second.state.isUsed ) {
						if ( _config.normalization == NORM_TRACE )
							*ofs << originTime.iso()
								 << " " << dit->second.state.fit
								 << " " << (dit->second.state.traceFit * traceScale)
								 << " " << traceScale
								 << std::endl;
						else
							*ofs << originTime.iso()
								 << " " << dit->second.state.fit
								 << " " << (dit->second.state.fitRaw * totalScale)
								 << " " << totalScale
								 << std::endl;
					}
					else {
						*ofs << originTime.iso()
							 << " " << dit->second.state.fit
							 << " 0 0 "
							 << std::endl;
					}
				}
			}

			if ( _config.detectorOffThreshold < 0 )
				--evt.processing.triggerCounter;
			else {
				if ( fit < _config.detectorThreshold )
					evt.processing.triggerCounter = 0;
			}

			if ( !evt.processing.triggerCounter ) {
				double timeOffset = 0.0;
				for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
					if ( !dit->second.dataState.isUsed ) continue;
					timeOffset += dit->second.dataState.timeLag;
				}

				timeOffset /= evt.processing.bestFitChannelCount;
				foundMatch(evt, timeOffset, stationCount, channelCount);
				evt.processing.reset();

				if ( evt.processing.triggerFitOutput.is_open() ) {
					// Close all trigger output
					evt.processing.triggerFitOutput.close();
					OutputMap::iterator oit;
					for ( oit = evt.processing.triggerFitChannelOutputs.begin();
					      oit != evt.processing.triggerFitChannelOutputs.end(); ++oit )
						oit->second.reset();
				}
			}
		}

		if ( _config.dumpFitFunction ) {
			evt.processing.fitOutput << evt.processing.loop << " " << fit
			                         << " " << originTime.iso() << " "
			                         << ((double)evt.processing.triggerCounter/(double)searchMaximumSteps)
			                         << std::endl;

			for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
				if ( !dit->second.dataState.trace ) continue;

				OutputMap::iterator oit = evt.processing.fitChannelOutputs.find(dit->first);
				ofstream_ptr ofs;
				if ( oit == evt.processing.fitChannelOutputs.end() ) {
					ofs = ofstream_ptr(new std::ofstream);
					evt.processing.fitChannelOutputs[dit->first] = ofs;
					ofs->open((evt.name + "-" + dit->first + ".fit").c_str());
				}
				else
					ofs = oit->second;

				if ( dit->second.state.isUsed ) {
					if ( _config.normalization == NORM_TRACE )
						*ofs << evt.processing.loop
						     << " " << dit->second.state.fit
						     << " " << (dit->second.state.traceFit * traceScale)
						     << " " << traceScale
						     << " " << originTime.iso() << std::endl;
					else
						*ofs << evt.processing.loop
						     << " " << dit->second.state.fit
						     << " " << (dit->second.state.fitRaw * totalScale)
						     << " " << totalScale
						     << " " << originTime.toString("%Y%m%d%H%M%S%6f") << std::endl;
				}
				else {
					*ofs << evt.processing.loop
					     << " " << dit->second.state.fit
					     << " 0 0 "
					     << originTime.toString("%Y%m%d%H%M%S%6f") << std::endl;
				}
			}
		}

		time += step;
		++evt.processing.loop;

		// Break if wrong step
		if ( step.seconds() < 0 || step.microseconds() < 0 ) break;
		if ( step.seconds() == 0 && step.microseconds() == 0 ) break;
	}

	// Remove temporary traces
	for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit )
		dit->second.dataState.trace = NULL;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
const Detector::ProcessingResult *
Detector::EventGroup::pending(const MasterEvent *evt,
                              const Core::Time &eventTime,
                              double margin) const {
	ResultQueue::const_iterator it;
	for ( it = queuedResults.begin(); it != queuedResults.end(); ++it ) {
		if ( it->event != evt ) continue;
		if ( it->time >= eventTime - Core::TimeSpan(margin) &&
		     it->time <= eventTime + Core::TimeSpan(margin) )
			return &(*it);
	}

	return NULL;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool Detector::EventGroup::popPending(ProcessingResult &res,
                                      const Core::Time &eventTime,
                                      double margin) {
	ResultQueue::iterator it;
	for ( it = queuedResults.begin(); it != queuedResults.end(); ) {
		if ( it->time >= eventTime - Core::TimeSpan(margin) &&
		     it->time <= eventTime + Core::TimeSpan(margin) ) {
			res = *it;
			queuedResults.erase(it);
			return true;
		}
	}

	return false;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::EventGroup::removePending(const Core::Time &eventTime, double margin) {
	ResultQueue::iterator it;
	for ( it = queuedResults.begin(); it != queuedResults.end(); ) {
		if ( it->time >= eventTime - Core::TimeSpan(margin) &&
			 it->time <= eventTime + Core::TimeSpan(margin) )
			it = queuedResults.erase(it);
		else
			++it;
	}
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::foundMatch(const MasterEvent &evt, double timeOffset,
                          int stationCount, int channelCount) {
	DataSet::const_iterator dit;

	std::cerr << "*=========================================================*" << std::endl;
	std::cerr << "  " << std::setprecision(2) << evt.processing.bestFit << " match with event " << evt.name << std::endl;
	std::cerr << "  " << evt.processing.bestFitStationCount << "/" << stationCount << " stations" << std::endl;
	std::cerr << "  " << evt.processing.bestFitChannelCount << "/" << channelCount << " channels" << std::endl;
	std::cerr << "*---------------------------------------------------------*" << std::endl;
	std::cerr << "  Time       : " << evt.processing.bestFittingTime.toString("%F %T.%6f") << std::endl;
	std::cerr << "  Hypocenter : " << evt.latitude << " " << evt.longitude << " " << evt.depth << "km" << std::endl;
	std::cerr << "  Magnitude  : " << evt.processing.bestFitMag << std::endl;
	std::cerr << "  Negative   : " << evt.doNotSendOrigins << std::endl;
	std::cerr << "*---------------------------------------------------------*" << std::endl;
	for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
		if ( !dit->second.dataState.isUsed ) continue;
		std::cerr << "  " << dit->first << "  " << std::setprecision(2) << dit->second.dataState.fit;
		std::cerr << "  " << std::setprecision(2) << dit->second.dataState.traceFit;
		//std::cerr << "  " << dit->second.dataState.timeLag << "s";
		std::cerr << std::endl;
	}
	std::cerr << "*=========================================================*" << std::endl;

	if ( _config.dumpWaveforms ) {
		Core::TimeSpan margin(_config.detectorMaxTimeLag);
		bool error = false;

		std::string path = _config.dumpWaveformsPath + evt.name;
		if ( !Util::pathExists(path) ) {
			if ( !Util::createPath(path) ) {
				SEISCOMP_ERROR("Unable to create waveform plot output directory: %s",
				               path.c_str());
				error = true;
			}
		}

		if ( !error ) {
			path += '/';

			for ( dit = evt.data.begin(); dit != evt.data.end(); ++dit ) {
				if ( !dit->second.dataState.isUsed ) continue;

				Core::TimeWindow eventTw(evt.processing.bestFittingTime+Core::TimeSpan(evt.waveParams.offset)-margin,
				                         evt.processing.bestFittingTime+Core::TimeSpan(evt.waveParams.offset)+margin + Core::TimeSpan(evt.waveParams.duration));

				Trace &trace = *dit->second.dataState.trace;
				const RTTrace &rt_data = *dit->second.dataState.source;

				// Extract data
				if ( !rt_data.extract(*dit->second.dataState.trace, eventTw) ) {
					SEISCOMP_WARNING("  %s: data extraction failed", dit->first.c_str());
					dit->second.state.isUsed = false;
					continue;
				}

				correct(*dit->second.dataState.trace, dit->second.dataState.gain);
				prepareTrace(*dit->second.dataState.trace, evt.waveParams, &dit->second.dataState.pgv);

				DoubleArray *data = DoubleArray::Cast(trace.data());
				const DoubleArray *master = DoubleArray::ConstCast(dit->second.trace.data());

				int cnt = std::min(data->size(), master->size());
				double *data_samples = data->typedData();
				const double *master_samples = master->typedData();

				dit->second.dataState.norm = norm(cnt, data_samples, dit->second.dataState.normSum, dit->second.dataState.lag);
				dit->second.state.norm = norm(cnt, master_samples, dit->second.state.normSum, -dit->second.dataState.lag);

				int d_ofs = dit->second.dataState.lag;
				int m_ofs = 0;
				if ( d_ofs < 0 ) {
					m_ofs = -d_ofs;
					d_ofs = 0;
					cnt -= m_ofs;
				}
				else if ( d_ofs > 0 )
					cnt -= d_ofs;

				double scaleD = 1.0 / dit->second.dataState.norm;
				double scaleM = 1.0 / dit->second.state.norm;

				Core::Time t = evt.processing.bestFittingTime;
				double dt = 1.0 / trace.samplingFrequency();

				std::ofstream out;
				out.open((path + evt.processing.bestFittingTime.toString("%Y%m%d%H%M%S%6f") + "-" + dit->second.trace.streamID() + ".plot").c_str());
				for ( int i = 0; i < cnt; ++i, t += dt )
					out << scaleD*data_samples[i+d_ofs]
					    << " " << scaleM*master_samples[i+m_ofs]
					    << " " << t.toString("%Y%m%d%H%M%S%6f")
					    << std::endl;
			}
		}
	}

	Core::Time originTime = evt.processing.bestFittingTime + Core::TimeSpan(timeOffset);
	bool publishEvent = true;
	const ProcessingResult *bestGroupEvent = NULL;
	EventGroup *group = NULL;

	// Check if other events of the same group have triggered already and if the
	// fit is the largest. If not, queue the result until all other events have
	// processed this time window as well or a timeout.
	if ( !evt.group.empty() ) {
		EventGroups::iterator it = _eventGroups.find(evt.group);
		if ( it == _eventGroups.end() ) {
			SEISCOMP_WARNING("Internal error: events for group '%s' not registered: bail out",
			                 evt.group.c_str());
			quit();
			return;
		}

		group = &it->second;

		MasterEvents::iterator eit;
		double bestGroupFit = evt.processing.bestFit;

		for ( eit = group->events.begin(); eit != group->events.end(); ++eit ) {
			MasterEvent *other = eit->get();
			if ( other == &evt ) continue;

			// Get buffered result for event 'other' and the current origin time
			const ProcessingResult *res = group->pending(other, originTime, _config.detectorTimeWindow);

			// No event pending for this time period?
			if ( res == NULL ) {
				// Time window has been processed already including the current
				// event -> no detection for this master event
				if ( other->processedTimeWindow.endTime() >= originTime+Core::TimeSpan(_config.detectorTimeWindow+other->waveParams.offset+other->waveParams.duration) ) {
					SEISCOMP_DEBUG("No events pending for %s and processed time window includes already "
					               "the current event", other->name.c_str());
					continue;
				}

				// Otherwise do not publish and wait for results
				publishEvent = false;
				SEISCOMP_INFO("Delay publication of event %s and wait for other events "
				              "in the group to be processed.", evt.name.c_str());
				break;
			}
			else {
				// Check for maximum fit
				if ( res->fit > bestGroupFit ) {
					bestGroupEvent = res;
					bestGroupFit = res->fit;
				}
			}
		}
	}

	ProcessingResult res;
	res.event = &evt;
	res.time = originTime;
	res.fit = evt.processing.bestFit;
	res.magnitude = evt.processing.bestFitMag;
	res.associatedChannelCount = channelCount;
	res.associatedStationCount = stationCount;
	res.usedChannelCount = evt.processing.bestFitChannelCount;
	res.usedStationCount = evt.processing.bestFitStationCount;
	res.usedChannels = evt.processing.bestFitChannels;
	res.timeout = _config.detectorPublicationTimeout;

	if ( !publishEvent ) {
		// Add to pending events
		group->queuedResults.push_back(res);
		return;
	}

	if ( bestGroupEvent != NULL ) {
		SEISCOMP_DEBUG("Found better queued event");
		publish(bestGroupEvent);
	}
	else
		publish(&res);

	// Remove all pending events for this origin time
	if ( group != NULL ) {
		group->removePending(originTime, _config.detectorTimeWindow);
		SEISCOMP_DEBUG("Events still pending in group %s: %d", evt.group.c_str(),
		               (int)group->queuedResults.size());
	}
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void Detector::publish(const ProcessingResult *res) {
	if ( res->event->doNotSendOrigins ) {
		SEISCOMP_INFO("Do not publish negative event %s at %s with fit of %d%%",
		              res->event->name.c_str(), res->time.toString("%F %T.%f").c_str(),
		              int((res->fit) * 100 + 0.5));
		return;
	}

	SEISCOMP_INFO("Publish event %s at %s with fit of %d%%",
	              res->event->name.c_str(), res->time.toString("%F %T.%f").c_str(),
	              int((res->fit) * 100 + 0.5));

	if ( !_config.outputEvents.empty() ) {
		std::ostream *os;
		std::ofstream ofs;

		if ( _config.outputEvents == "-" )
			os = &std::cout;
		else {
			ofs.open(_config.outputEvents.c_str(), std::ios_base::out | std::ios_base::app);
			os = &ofs;
		}

		// yyyy mm dd HH MM SS.FFF Lat Long Mag city CF
		*os << res->time.toString("%Y %m %d %H %M %S.%3f")
		    << "  " << std::fixed << std::setprecision(3) << res->event->latitude
		    << "  " << std::fixed << std::setprecision(3) << res->event->longitude
		    << "  " << std::fixed << std::setprecision(1) << res->magnitude
		    << "  " << (res->event->place.empty()?"-":res->event->place)
		    << "  " << std::fixed << std::setprecision(2) << res->fit
		    << "  " << res->usedChannels.size() << "  (";

		ResultItems::const_iterator rit;
		bool first = true;
		for ( rit = res->usedChannels.begin(); rit != res->usedChannels.end(); ++rit ) {
			const ResultItem &item = *rit;
			if ( !first ) *os << ", ";
			*os << item.data->trace.streamID() << ":" << std::setprecision(2) << item.fit;
			first = false;
		}

		*os << ")" << std::endl;
	}

	//if ( _config.test ) return;

	SEISCOMP_DEBUG("Sending origin");

	Core::Time now = Core::Time::GMT();
	DataModel::CreationInfo ci;
	ci.setAgencyID(agencyID());
	ci.setAuthor(author());
	ci.setCreationTime(now);

	DataModel::MagnitudePtr mag = DataModel::Magnitude::Create();
	mag->setCreationInfo(ci);
	mag->setType("MAGS");
	mag->setMagnitude(DataModel::RealQuantity(res->magnitude));
	mag->setStationCount(res->usedStationCount);

	DataModel::OriginPtr org = DataModel::Origin::Create();
	org->setCreationInfo(ci);
	org->setLatitude(DataModel::RealQuantity(res->event->latitude));
	org->setLongitude(DataModel::RealQuantity(res->event->longitude));
	org->setDepth(DataModel::RealQuantity(res->event->depth));
	org->setTime(DataModel::TimeQuantity(res->time));
	org->setEpicenterFixed(true);
	org->setEvaluationMode(DataModel::EvaluationMode(DataModel::AUTOMATIC));

	std::vector<double> azis;
	std::vector<double> dists;
	DataSet::const_iterator dit;

	for ( dit = res->event->data.begin(); dit != res->event->data.end(); ++dit ) {
		double az, baz, dist;
		Math::Geo::delazi(res->event->latitude, res->event->longitude,
		                  dit->second.stationData->model->latitude(),
		                  dit->second.stationData->model->longitude(),
		                  &dist, &az, &baz);

		dists.push_back(dist);
		azis.push_back(az);
	}

	std::sort(azis.begin(), azis.end());
	std::sort(dists.begin(), dists.end());

	DataModel::OriginQuality oq;

	if ( azis.size() > 2 ) {
		double azGap = 0;
		for ( size_t i = 0; i < azis.size()-1; ++i )
			azGap = (azis[i+1]-azis[i]) > azGap ? (azis[i+1]-azis[i]) : azGap;

		oq.setAzimuthalGap(azGap);
		mag->setAzimuthalGap(azGap);
	}

	if ( !dists.empty() ) {
		oq.setMinimumDistance(dists.front());
		oq.setMaximumDistance(dists.back());
		oq.setMedianDistance(dists[dists.size()/2]);
	}

	oq.setStandardError(1.0-res->fit);
	oq.setAssociatedStationCount(res->associatedStationCount);
	oq.setUsedStationCount(res->usedStationCount);
	oq.setAssociatedPhaseCount(res->associatedChannelCount);
	oq.setUsedPhaseCount(res->usedChannelCount);

	org->setQuality(oq);
	org->setMethodID("MAGS");

	DataModel::Notifier::Enable();
	DataModel::EventParameters ep;

	// Create the objects top-down
	ep.add(org.get());
	org->add(mag.get());

	for ( size_t i = 0; i < res->event->picks.size(); ++i ) {
		DataModel::PickPtr p = DataModel::Pick::Create();
		p->setCreationInfo(ci);

		// Apply time offset according to new detection
		Core::TimeSpan ofs = res->event->picks[i].pick->time().value() - res->event->time;
		p->setTime(res->time + ofs);
		p->setWaveformID(res->event->picks[i].pick->waveformID());
		p->setEvaluationMode(DataModel::EvaluationMode(DataModel::AUTOMATIC));

		try { p->setPhaseHint(res->event->picks[i].pick->phaseHint()); }
		catch ( ... ) {}

		// Add pick to output
		ep.add(p.get());

		// Create arrival
		DataModel::ArrivalPtr ar = new DataModel::Arrival;
		ar->setPickID(p->publicID());
		ar->setPhase(res->event->picks[i].phase);
		if ( res->event->picks[i].weight )
			ar->setWeight(*res->event->picks[i].weight);

		// Attach arrival to origin
		org->add(ar.get());
	}

	ResultItems::const_iterator rit;
	for ( rit = res->usedChannels.begin(); rit != res->usedChannels.end(); ++rit ) {
		const ResultItem &item = *rit;
		// Create station magnitude
		DataModel::StationMagnitudePtr stamag = DataModel::StationMagnitude::Create();
		stamag->setCreationInfo(ci);
		stamag->setMagnitude(item.magnitude);
		stamag->setOriginID(org->publicID());
		stamag->setType("MAGS");
		stamag->setWaveformID(DataModel::WaveformStreamID(item.data->trace.networkCode(),
		                                                  item.data->trace.stationCode(),
		                                                  item.data->trace.locationCode(),
		                                                  item.data->trace.channelCode(),
		                                                  ""));

		// Add reference to network magnitude
		DataModel::StationMagnitudeContributionPtr ref = new DataModel::StationMagnitudeContribution;
		ref->setStationMagnitudeID(stamag->publicID());
		ref->setWeight(1.0);
		mag->add(ref.get());
		org->add(stamag.get());
	}

	DataModel::NotifierMessagePtr msg = DataModel::Notifier::GetMessage();

	DataModel::Notifier::Disable();

	logObject(_logOrigins, now);

	if ( connection() && msg )
		connection()->send(msg.get());
	/*
	else {
		IO::XMLArchive ar;
		DataModel::EventParameters *tmp = &ep;
		ar.create("-");
		ar.setFormattedOutput(true);
		ar << tmp;
		ar.close();
	}
	*/
}
