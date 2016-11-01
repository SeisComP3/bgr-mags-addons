/***************************************************************************
 *   Copyright (C) by GFZ Potsdam                                          *
 *                                                                         *
 *   You can redistribute and/or modify this program under the             *
 *   terms of the SeisComP Public License.                                 *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   SeisComP Public License for more details.                             *
 ***************************************************************************/

#include <iostream>
#include <fstream>

#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/io/records/sacrecord.h>
#include <seiscomp3/core/recordsequence.h>

#include "filter.h"

using namespace Seiscomp;


int main(void)
{
	IO::RecordStreamPtr rs = IO::RecordStream::Create("file");
	rs->setSource("-");
	IO::RecordInput inp(rs.get(), Array::DOUBLE, Record::DATA_ONLY);

	RecordPtr rec;
	RingBuffer seq(Core::TimeSpan(6,0));

	while ( rec = inp.next() )
		seq.feed(rec.get());

	RecordPtr r = seq.continuousRecord<double>();
	const DoubleArray *ar = DoubleArray::ConstCast(r->data());

	std::vector<double> in, out;

	in.resize(ar->size());

	double mean = ar->mean();

	for ( int i = 0; i < ar->size(); ++i )
		in[i] = (*ar)[i] - mean;

	std::vector<Math::Complex> spec;

	out = in;

	Math::fft(spec, in);
	bandpass(spec, 1.0/r->samplingFrequency(), 3, 10, 40);
	hilbert(spec);
	Math::ifft(out, spec);
	magnitude(&in[0], in.size(), &out[0]);
	//envelope(&in[0], in.size(), &out[0]);

	for ( int i = 0; i < ar->size(); ++i )
		std::cout << in[i] << "\t" << out[i] << "\t" << std::endl;

	IO::SACRecord sac(r->networkCode(), r->stationCode(), r->locationCode(),
	                  r->channelCode(), r->startTime(), r->samplingFrequency());

	std::ofstream ofs;

	sac.setData(in.size(), &in[0], Array::DOUBLE);
	ofs.open("in.sac");
	sac.write(ofs);
	ofs.close();

	sac.setLocationCode("00");
	sac.setData(out.size(), &out[0], Array::DOUBLE);
	ofs.open("env-jan.sac");
	sac.write(ofs);
	ofs.close();
}
