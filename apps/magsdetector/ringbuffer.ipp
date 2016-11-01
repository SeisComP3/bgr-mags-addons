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


#include <string.h>


template <typename T>
CircularBuffer<T>::CircularBuffer(size_t capacity)
: _capacity(0), _start(0), _end(0), _data(NULL)  {
	setCapacity(capacity);
}


template <typename T>
CircularBuffer<T>::~CircularBuffer() {
	if ( _data != NULL ) delete _data;
}


template <typename T>
void CircularBuffer<T>::setCapacity(size_t n) {
	// No change in capacity
	if ( _capacity == n ) return;

	// Delete existing data block
	if ( _data != NULL ) {
		delete _data;
		_data = NULL;
	}

	_capacity = n;
	if ( _capacity > 0 )
		_data = new T[_capacity+1];

	// Reset all index pointers
	_start = _end = 0;
}


template <typename T>
size_t CircularBuffer<T>::capacity() const {
	return _capacity;
}


template <typename T>
size_t CircularBuffer<T>::size() const {
	if ( _end < _start ) return _capacity-_start+_end+1;
	return _end-_start;
}


template <typename T>
bool CircularBuffer<T>::empty() const {
	return _start == _end;
}


template <typename T>
void CircularBuffer<T>::append(const T *d, size_t n) {
	if ( n >= _capacity ) {
		d += n-_capacity;
		n = _capacity;
		_start = 0;
		_end = _capacity;
		memcpy(_data, d, _capacity*sizeof(T));
		return;
	}

	size_t new_end = _end + n;
	bool move_front = size() + n > _capacity;
	if ( new_end > _capacity ) {
		new_end -= _capacity;
		--new_end;
	}

	// Simple case
	if ( new_end > _end )
		memcpy(_data + _end, d, n*sizeof(T));
	else {
		// Split up blocks
		size_t remain = _capacity - _end + 1;
		memcpy(_data + _end, d, remain*sizeof(T));
		n -= remain;
		d += remain;
		memcpy(_data, d, n*sizeof(T));
	}

	_end = new_end;
	if ( move_front ) _start = _end+1;
}


template <typename T>
void CircularBuffer<T>::append(const T *d, size_t n, T scale) {
	if ( n >= _capacity ) {
		d += n-_capacity;
		n = _capacity;
		_start = 0;
		_end = _capacity;
		memcpy(_data, d, _capacity*sizeof(T));
		for ( size_t i = 0; i < _capacity; ++i ) _data[i] *= scale;
		return;
	}

	size_t new_end = _end + n;
	bool move_front = size() + n > _capacity;
	if ( new_end > _capacity ) {
		new_end -= _capacity;
		--new_end;
	}

	// Simple case
	if ( new_end > _end ) {
		T *block = _data + _end;
		memcpy(block, d, n*sizeof(T));
		for ( size_t i = 0; i < n; ++i ) block[i] *= scale;
	}
	else {
		// Split up blocks
		size_t remain = _capacity - _end + 1;
		T *block = _data + _end;
		memcpy(block, d, remain*sizeof(T));
		for ( size_t i = 0; i < remain; ++i ) block[i] *= scale;
		n -= remain;
		d += remain;
		memcpy(_data, d, n*sizeof(T));
		for ( size_t i = 0; i < n; ++i ) _data[i] *= scale;
	}

	_end = new_end;
	if ( move_front ) _start = _end+1;
}


template <typename T>
size_t CircularBuffer<T>::copy(T *target, size_t ofs, size_t n) const {
	size_t copied = 0;
	size_t size_ = size();

	if ( ofs >= size_ ) return copied;

	size_ -= ofs;

	size_t s = _start + ofs;
	if ( s > _capacity ) {
		s -= _capacity;
		--s;
	}

	// Simple case
	if ( _end > s ) {
		size_t to_copy = std::min(n, _end-s);
		memcpy(target, _data+s, to_copy*sizeof(T));
		copied += to_copy;
	}
	else {
		// Split up blocks
		size_t to_copy = std::min(n, _capacity - s + 1);
		memcpy(target, _data+s, to_copy*sizeof(T));
		copied += to_copy;
		n -= to_copy;

		if ( n > 0 ) {
			target += to_copy;
			to_copy = std::min(n, _end);
			memcpy(target, _data, to_copy*sizeof(T));
			copied += to_copy;
		}
	}

	return copied;
}


template <typename T>
void CircularBuffer<T>::clear() {
	_start = _end = 0;
}


template <typename T>
void CircularBuffer<T>::dump(std::ostream &os) {
	os << "Capacity: " << capacity() << std::endl;
	os << "Size    : " << size() << std::endl;
	os << "get_ptr : " << _start << std::endl;
	os << "put_ptr : " << _end << std::endl;

	size_t i = _start;
	while ( i != _end ) {
		if ( i != _start ) os << ",";
		os << _data[i];
		++i; if ( i > _capacity ) i = 0;
	}
	os << std::endl;
}
