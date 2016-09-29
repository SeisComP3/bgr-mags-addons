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


#ifndef __MAGS_DETECTOR_H__
#define __MAGS_DETECTOR_H__


#include "app.h"
#include "ringbuffer.h"

#include <seiscomp3/core/datetime.h>
#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/recordsequence.h>
#include <seiscomp3/core/typedarray.h>
#include <seiscomp3/datamodel/pick.h>

#include <fstream>
#include <stdint.h>


namespace Seiscomp {


class Detector : public Application {
	public:
		DEFINE_SMARTPOINTER(Trace);
		struct Trace : GenericRecord {
			int    signalOffset;
			int    signalCnt;
			int    noiseOffset;
			int    noiseCnt;
			double mean;

			const double *signal() const { return DoubleArray::ConstCast(data())->typedData() + signalOffset; }
			const double *noise() const { return DoubleArray::ConstCast(data())->typedData() + noiseOffset; }

			int signalCount() const { return signalCnt; }
			int noiseCount() const { return noiseCnt; }
		};

		DEFINE_SMARTPOINTER(RTTrace);
		class RTTrace : public Core::BaseObject {
			public:
				RTTrace() : _samplingFrequency(0) {}

				void setSamplingFrequency(double sf);
				double samplingFrequency() const { return _samplingFrequency; }

				void setGain(double gain) { _gain = gain; }

				void setTimeSpan(const Core::TimeSpan &);

				void reset();
				void feed(Record *rec);

				const Core::Time &startTime() const {
					return _startTime;
				}

				const Core::Time &endTime() const {
					return _endTime;
				}

				size_t sampleCount() const {
					return _data.size();
				}

				bool extract(Trace &trace, const Core::TimeWindow &tw) const;

			private:
				CircularBuffer<double>   _data;
				Core::Time               _startTime;
				Core::Time               _endTime;
				double                   _samplingFrequency;
				double                   _gain;

				Core::TimeSpan           _maxAllowedGap;
				Core::TimeSpan           _maxAllowedOverlap;
		};

		struct State {
			bool        isUsed;
			double      pgv;
			double      norm;
			double      normSum;
			double      fit;
			double      fitRaw;
			double      timeLag;
			int         lag;
			double      traceFit;
		};

		struct DataState : State {
			RTTrace    *source;
			TracePtr    trace;
			double      gain;
		};

		DEFINE_SMARTPOINTER(StationData);
		struct StationData : public Core::BaseObject {
			StationData(DataModel::SensorLocation *l) : model(l) {}

			DataModel::SensorLocation *model;
			int                        channelsUsed;
		};

		struct Data {
			Data() : seq(NULL) {}

			void dispose() {
				if ( seq != NULL ) {
					delete seq;
					seq = NULL;
				}
			}

			Trace                     trace;
			RecordSequence           *seq;

			mutable State             state;
			mutable DataState         dataState;

			mutable StationDataPtr    stationData;
		};


		typedef std::map<std::string, Data> DataSet;
		typedef std::map<std::string, StationDataPtr> StationSet;
		typedef std::map<std::string, RTTracePtr> Waveforms;

		typedef std::pair<double,double> FilterFreqs;
		typedef std::vector<std::string> StringList;

		typedef std::list<DataModel::WaveformStreamID> WaveformIDList;
		typedef std::list<DataModel::Stream*> StreamList;

		enum Normalization {
			NORM_TRACE,
			NORM_TOTAL
		};

		struct WaveParams {
			int           filterOrder;
			FilterFreqs   filterFreqs;
			bool          filterBandStop;

			bool          envelopeEnable;
			int           envelopeSamplingFrequency;
			bool          envelopeDownsampleAverage;
			double        envelopeHiFreq;
			bool          envelopeAcausal;

			double        noiseBegin;
			double        noiseEnd;
			double        noise2Begin;
			double        noise2End;
			double        signalBegin;
			double        signalEnd;

			double        offset;
			double        duration;

			bool          processingLog;
			bool          processingAcausal;
		};

		struct Config {
			Config();

			bool          offline;
			bool          test;
			bool          dumpFitFunction;
			int           dumpFitChannels;

			StringList    channels;
			StringList    events;

			bool          dumpWaveforms;
			bool          dumpTriggerWaveforms;
			bool          dumpTriggerFit;
			std::string   dumpWaveformsPath;

			int           bufferSize;
			int           processingInterval;

			WaveParams    waveParams;
			double        maximumLatency;

			std::string   normalizationString;
			Normalization normalization;

			double        detectorThreshold;
			double        detectorOffThreshold;
			double        detectorChannelThreshold;
			double        detectorTimeWindow;
			double        detectorMaxTimeLag;
			int           detectorMinStationRatio;
			int           detectorMinChannelRatio;
			double        detectorProcessingTimeWindow;
			int           detectorPublicationTimeout;
			int           detectorMaximumStepFrequency;

			std::string   outputEvents;
			bool          showProgress;
		};

		typedef boost::shared_ptr<std::ofstream> ofstream_ptr;
		typedef std::map<std::string, ofstream_ptr> OutputMap;

		struct ResultItem {
			ResultItem(const Data *d, double f, double m)
			: data(d), fit(f), magnitude(m) {}
			const Data  *data;
			double       fit;
			double       magnitude;
		};

		typedef std::vector<ResultItem> ResultItems;

		struct ProcessingState {
			ProcessingState() : loop(0) {}
			uint64_t         loop;
			Core::Time       bestFittingTime;
			double           bestFit;
			double           bestFitMag;
			int              bestFitStationCount;
			int              bestFitChannelCount;
			ResultItems      bestFitChannels;

			int              triggerCounter;

			std::ofstream    fitOutput;
			OutputMap        fitChannelOutputs;

			std::ofstream    triggerFitOutput;
			OutputMap        triggerFitChannelOutputs;

			void reset() {
				bestFittingTime = Core::Time();
				bestFit = 0;
				bestFitStationCount = -1;
				bestFitChannelCount = -1;
				bestFitChannels.clear();
				triggerCounter = -1;
			}
		};

		struct EventArrival {
			DataModel::PickPtr pick;
			DataModel::Phase   phase;
			OPT(double)        weight;
		};

		typedef std::vector<EventArrival> PickList;

		DEFINE_SMARTPOINTER(MasterEvent);
		struct MasterEvent : Core::BaseObject {
			std::string             name;
			std::string             place;
			std::string             group;
			Core::Time              time;
			double                  magnitude;
			OPT(double)             deltaM;
			double                  latitude;
			double                  longitude;
			double                  depth;
			std::string             dataURI;
			DataSet                 data;
			StationSet              stations;
			bool                    processData;
			bool                    doNotSendOrigins;

			WaveParams              waveParams;
			PickList                picks;

			Core::TimeWindow        processedTimeWindow;

			mutable ProcessingState processing;

			bool read(const Seiscomp::Config::Config &cfg,
			          DataModel::DatabaseQuery *db,
			          const std::string &prefix);
			bool readData(const Config &cfg, const StreamList &streams);
			bool prepareData(const Config &);
			void dumpConfig(std::ostream &os) const;
		};

		typedef std::list<MasterEventPtr> MasterEvents;

		struct ProcessingResult {
			const MasterEvent *event;
			Core::Time         time;
			double             fit;
			double             magnitude;
			int                associatedStationCount;
			int                associatedChannelCount;
			int                usedStationCount;
			int                usedChannelCount;
			ResultItems        usedChannels;
			int                timeout;
		};

		typedef std::list<ProcessingResult> ResultQueue;

		struct EventGroup {
			MasterEvents events;
			ResultQueue  queuedResults;

			const ProcessingResult *pending(const MasterEvent *evt, const Core::Time &eventTime,
			                                double margin) const;
			bool popPending(ProcessingResult &res, const Core::Time &eventTime,
			                double margin);
			void removePending(const Core::Time &eventTime, double margin);
		};

		typedef std::map<std::string, EventGroup> EventGroups;

	public:
		Detector(int argc, char **argv);
		~Detector();

	protected:
		void printVersion();
		bool validateParameters();

		bool init();
		bool run();
		void done();

		void handleRecord(Record *rec);
		void handleTimeout();


	private:
		void resolveStreams(StreamList &streams, const Core::Time &);
		void dumpWaveformState();
		void dumpMasterEvents();
		void dumpWaveforms(const MasterEvent &evt, const Core::Time &,
		                   double fit);
		void processWaveforms();

		static bool prepareTrace(Trace &trace, const WaveParams &params,
		                         double *pgv);

		static bool prepareSignalWindow(Trace &trace, const WaveParams &params);

		void process(const MasterEvent &evt, Core::TimeWindow tw);
		void process(double minTimeWindow, double latency);

		void foundMatch(const MasterEvent &evt, double timeOffset,
		                int stationCount, int channelCount);

		void publish(const ProcessingResult *res);
		void flushResults(bool force);


	private:
		Config           _config;
		MasterEvents     _masterEvents;
		EventGroups      _eventGroups;
		WaveformIDList   _channels;
		Waveforms        _waveforms;

		Core::TimeWindow _lastProcessedTimeWindow;
		int              _lastProcessedChannelCount;
		int              _currentProcessingTicker;

		double           _minimumRequiredWindowLength;
		bool             _recordsReceived;

		ObjectLog       *_logOrigins;
};


}


#endif
