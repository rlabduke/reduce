// Name: VectIter.h
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: Vector iterator interface
//          (includes non-const vector iterator)
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#ifndef VECTITER_H
#define VECTITER_H 1

#include "Iter.h"
#include "Vector.h"

#ifndef BOOLPREDEFINED
typedef int bool;
#endif

template <class T>
class VectIter : public Iterator<T> {
public:
   VectIter(const Vector<T>& v) : _elem(0), _vec(v) {}
   VectIter(const VectIter<T>& vi) : _elem(vi._elem), _vec(vi._vec) {}

   void sync(const VectIter<T>& vi) { _elem = vi._elem; }

   // check current entry
   bool valid() const { return (_elem >= 0)
                            && (_elem < _vec._num_elem); }
   operator bool() const { return valid(); }

   int size() const { return _vec.size(); }

   const T& operator *() const{      // data access
      assert(valid());
      return _vec[_elem];
   }

   void reset() { _elem = 0; } // go to first entry
   void begin() { _elem = 0; }
   void end()   { _elem = _vec.size()-1; } // go to last entry

   bool next(T& e) {                   // return the entry and advance
      bool rc = valid();
      if (rc) { e = _vec[_elem]; _elem++; }
      return rc;
   }
   bool prev(T& e) {                   // return the entry and retreat
      bool rc = valid();
      if (rc) { e = _vec[_elem]; _elem--; }
      return rc;
   }
                             // retreat to prev entry                             // advance to next entry
   bool operator++()    { _elem++; return valid(); }
   bool operator++(int) { _elem++; return valid(); }

   bool operator--()    { _elem--; return valid(); }
   bool operator--(int) { _elem--; return valid(); }
protected:
   VectIter<T>& operator=(const VectIter<T>&); // can't assign

   int   _elem;            // current element
   const Vector<T>& _vec;  // the Vector
};

template <class T>
class NonConstVectIter : public VectIter<T> {
public:
   NonConstVectIter(Vector<T>& v) : VectIter<T>(v), _NCvec(v) {}
   NonConstVectIter(const NonConstVectIter<T>& vi)
			      : VectIter<T>(vi), _NCvec(vi._NCvec) {}

   T& data() const{                  // non-const data access
      assert(valid());
      return _NCvec[_elem];
   }

   void update(const T& e) const{    // data modication
      assert(valid());
      _NCvec[_elem] = e;
   }
protected:
   NonConstVectIter<T>& operator=(const NonConstVectIter<T>&); // can't assign

   Vector<T>& _NCvec;  // non-const reference to the Vector
};
#endif
