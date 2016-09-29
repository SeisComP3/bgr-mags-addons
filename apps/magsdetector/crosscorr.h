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


#ifndef __MAGS_CROSSCORR_H__
#define __MAGS_CROSSCORR_H__


#include <math.h>
#include <climits>


template <typename T, class NORM>
void goodnessOfFit(const T *x, int n, const T *y, NORM &norm);

template <typename T, class NORM>
int crosscorr(const T *x, int n, const T *y, int maxDelay, NORM &norm);


#include "crosscorr.ipp"

#endif
