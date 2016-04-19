// Name: Vector.C
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: Vector implementation
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#include "Vector.h"
#include "utility.h" // for swap(x, y)

// allocates memory for elements of Vector<T>
template <class T>
void Vector<T>::buildArray(int n) {
   assert(n >= 0);

   _array = new T[n];
   _num_elem  = n;
   _array_len = n;
}

template <class T>
Vector<T>::Vector(int n) {
   buildArray(n);
}

template <class T>
Vector<T>::Vector(int sz, const T& defaultVal) {
   buildArray(sz);
   for (int i=0; i < sz; i++) { _array[i] = defaultVal; }
}

// copy constructor
template <class T>
Vector<T>::Vector(const Vector& v) {
   buildArray(v.size());
   for (int i=0; i < v.size(); i++) {
      _array[i] = v._array[i];
   }
}

// public member functions:

// returns element value (non-const object)
template <class T>
T& Vector<T>::data(int index) {
   assert(index >= 0 && index < _num_elem);
   return _array[index];
}

// returns element value (non-const object)
template <class T>
T& Vector<T>::operator[](int index) {
   assert(index >= 0 && index < _num_elem);
   return _array[index];
}

// assign a new value
template <class T>
void Vector<T>::update(int index, const T& val) {
   assert (index >= 0 && index < _num_elem);
   _array[index] = val;
}

// lazy resize() -- only grows array unless new_n is negative
// change number of elements, maintaining existing data in range
template <class T>
void Vector<T>::resize(int new_n) {
   if (new_n < 0 || new_n > _array_len) {
      if (new_n < 0) { new_n = -new_n; }
      int cpycnt = (new_n < _num_elem) ? new_n : _num_elem;

      T* temp_array = new T[new_n];

      for (int i=0; i < cpycnt; i++) { // copy over old data
	 temp_array[i] = _array[i];
      }
      delete [] _array;
      _array     = temp_array;
      _array_len = new_n;
   }
   else if (new_n > _num_elem) {
      for (int i=_num_elem; i < new_n; i++) {
	 _array[i] = T();  // re-initialize
      }
   }
   _num_elem = new_n;
}

// heapsort (working indexes are 1 based)
template <class T>
void Vector<T>::sort(int len, int offset) {
   if (offset < 0 || offset >= _num_elem) { offset = 0; }
   if (len < 0 || len+offset > _num_elem) { len = _num_elem - offset; }

   int i;                                 // build heap
   for(i = len >> 1; i >= 1; i--) { // starting in the center
      siftdown(i, len, offset);
   }
   for(i = len; i >= 2; i--) { // put max values at the end
      swap(_array[offset], _array[i-1+offset]);
      siftdown(1, i-1, offset);
   }
}

// siftdown is a heap building utility func.
template <class T>
void Vector<T>::siftdown(int l, int u, int offset) {
//    (working indexes are 1 based)
//  pre-condition: maximal heap(l+1, u)
// post-condition: maximal heap(l,   u)

   int i = l;
   for(;;) { // maxheap(l, u) except between i and its children
      int c = i << 1; // double
      if (c > u) { break; }
      if (c+1 <= u && _array[c-1+offset] < _array[c+offset]) { c++; }
      if (! (_array[i-1+offset] < _array[c-1+offset])) { break; }
      swap(_array[c-1+offset], _array[i-1+offset]);
      i = c;
   }
}

// reverse the order of items
template <class T>
void Vector<T>::rev(int len, int offset) {
   if (offset < 0 || offset >= _num_elem) { offset = 0; }
   if (len < 0 || len+offset > _num_elem) { len = _num_elem - offset; }

   for(int i = len >> 1; i >= 1; i--) {
      swap(_array[i-1+offset], _array[len-i+offset]);
   }
}

// assignment operator
template <class T>
Vector<T>& Vector<T>::operator = (const Vector& v) {

   if (this != &v) {
      resize(v.size());

      for (int i=0; i < v.size(); i++) {
	 _array[i] = v._array[i];
      }
   }
   return *this ;
}

// subscript operator [] (const object)
template <class T>
const T& Vector<T>::operator [] (int index) const {
   assert (index >= 0 && index < _num_elem);
   return _array[index];
}

// concatenate
template <class T>
Vector<T>& Vector<T>::operator += (const Vector& v) {
   int offset = size();
   int n      = v.size();
   resize(size() + v.size());

   for (int i=0; i < n; i++) {
      update(offset+i, v._array[i]);
   }
   return *this;
}

// output operator
template <class T>
ostream& operator << (ostream& s, const Vector<T>& v) {
   s << "[" ;
   for (int i=0; i < v._num_elem; i++) {
      s << v[i] << ((i != v._num_elem - 1) ? ", " : "");
   }
   s << "]";
   return s ;
}

