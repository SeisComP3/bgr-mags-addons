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


#ifndef __MAGS_CIRCULARBUFFER_H__
#define __MAGS_CIRCULARBUFFER_H__


#include <ostream>


template <typename T>
class CircularBuffer {
	public:
		CircularBuffer(size_t capacity = 0);
		~CircularBuffer();

		// Sets the current capacity. Causes clear.
		void setCapacity(size_t n);

		size_t capacity() const;
		size_t size() const;
		bool empty() const;

		void clear();

		void append(const T *d, size_t n);
		void append(const T *d, size_t n, T scale);

		size_t copy(T *target, size_t ofs, size_t n) const;

		void dump(std::ostream &os);

		size_t _start_index() const { return _start; }
		size_t _end_index() const { return _end; }

	private:
		size_t  _capacity;
		size_t  _start;
		size_t  _end;
		T      *_data;
};


#include "ringbuffer.ipp"


#endif
