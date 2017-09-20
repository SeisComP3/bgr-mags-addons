<pre>
/***************************************************************************
* Copyright (C) by                                                         *
* Bundesanstalt fuer Geowissenschaften und Rohstoffe (BGR)                 *
* Geozentrum Hannover                                                      *
* Stilleweg 2                                                              *
* 30655 Hannover                                                           *
* Germany                                                                  *
*                                                                          *
* You can redistribute and/or modify this program under the                *
* terms of the SeisComP Public License.                                    *
*                                                                          *
* This program is distributed in the hope that it will be useful,          *
* but WITHOUT ANY WARRANTY; without even the implied warranty of           * 
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             *
* SeisComP Public License for more details.                                *
*                                                                          *
* This software was developed as part of the research project MAGS         *
* sposored by the Federal Ministry for Economic Affairs and Energy         *
* (FKZ: 0325191A)                                                          *
*                                                                          *
* When citing this software, please refer to the following article:        *
* M. Vasterling, U. Wegler, A. Bruestle, M. Bischoff, J. Becker, 2016.     *
* Real-time envelope cross-correlation detector: application to induced    *
* seismicity in the Insheim and Landau deep geothermal reservoirs.         *
* Journal of Seismology,  p1-16, issn: 1573-157X,                          *
* doi: 10.1007/s10950-016-9597-1                                           *
*                                                                          *
* Developed by gempa GmbH                                                  *
***************************************************************************/
</pre>
**Description**

The *magsdetector* detects events by finding similarities in defined master events based on cross-correlation of waveforms or envelopes. An event (origin) is declared if the waveform fit is above a configured threshold. Location and depth are taken from the master events whereas time is taken from the time window currently processed and magnitude is computed from the ratio of the master event PGV to the continuous PGV.

It can be ran in realtime or offline.

**Compile**

In order to use this module the sources have to be compiled to an executable. Merge them into the Seiscomp3 sources and compile Seiscomp3 as usual.
<pre>
# merge bgr-mags-addons and Seiscomp3
git clone https://github.com/SeisComP3/seiscomp3.git sc3-src
cd sc3-src
git submodule add -f https://github.com/SeisComP3/bgr-mags-addons.git src/bgr
</pre>
For compiling Seiscomp3, please refer to https://github.com/SeisComP3/seiscomp3#compiling
