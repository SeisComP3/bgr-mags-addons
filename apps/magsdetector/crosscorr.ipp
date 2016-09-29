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


template <typename T, class NORM>
int crosscorr(const T *x, int n, const T *y, int maxDelay, NORM &norm) {
	int i, j, delay, max_delay;
	T max_r;

	max_delay = 0;
	max_r = -std::numeric_limits<T>::max();

	/* Calculate the correlation series */
	for ( delay = -maxDelay; delay <= maxDelay; ++delay ) {
		NORM norm_(norm.scaleSyn, norm.scaleData);

		for ( i = 0; i < n; ++i ) {
			j = i + delay;
			if ( j < 0 || j >= n ) continue;
			norm_.update(x[i], y[j]);
		}

		T r = norm_();
		if ( r > max_r ) {
			max_delay = delay;
			max_r = r;
			norm = norm_;
		}
	}

	return max_delay;
}


template <typename T, class NORM>
void goodnessOfFit(const T *x, int n, const T *y, NORM &norm) {
	for ( int i = 0; i < n; ++i )
		norm.update(x[i], y[i]);
}
