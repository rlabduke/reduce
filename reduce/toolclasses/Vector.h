/// Name: Vector.h
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: Vector interface
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#ifndef VECTOR_H
#define VECTOR_H 1

#ifdef OLD_STD_HDRS
#include <assert.h>
#else
#include <cassert>
#endif
#include <iostream>

#ifdef BRACKETOPERPARMS
#define PARM_TYPE_T <T>
#else
#define PARM_TYPE_T
#endif

template <class T> class VectIter;
template <class T> class NonConstVectIter;

// indexes run from 0 to v.size()-1

template <class T>
class Vector {
   friend class VectIter<T>;
   friend class NonConstVectIter<T>;
public:
   Vector() { buildArray(0); }  // default size is 0
   Vector(int sz);
   Vector(int sz, const T& defaultVal);
   Vector(const Vector& v);
  ~Vector() { delete [] _array; _array=NULL; _num_elem=_array_len=0; }

   int size() const { return _num_elem; }

   const T& operator [](int index) const;
         T& operator [](int index); // vec[] can be an l-value
         T& data(int index);

   void update(int i, const T& val); // set the ith element to val

   // operator* returns the first value
   const T& operator*() { assert(_num_elem > 0); return _array[0]; }

   void resize(int new_n); // change the size of the array

   void sort(int len=-1, int offset=0); // heap sort low to high
   void rev (int len=-1, int offset=0); // reverse the order of values

   Vector& operator = (const Vector& v); // assignment

   Vector& operator += (const Vector& v); // concatenate

   friend std::ostream& operator << PARM_TYPE_T(std::ostream& s, const Vector<T>& v);
private:

   T*  _array;
   int _num_elem;
   int _array_len; // actual number of elements the array can hold 
		
   void buildArray(int n);
   void siftdown(int l, int u, int offset=0); // heap building utility func.
};

template <class T> std::ostream& operator << (std::ostream& s, const Vector<T>& v);

#ifdef INCTEMPLATEDEFNS
#include "Vector.cpp"
#endif

#endif
