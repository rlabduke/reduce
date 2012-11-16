// Name: Iterator.h
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: abstract iteration class
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#ifndef ITER_H 
#define ITER_H 1

#ifndef BOOLPREDEFINED
typedef int bool;
#endif

template <class T>
class Iterator {
public:
   virtual bool     valid()      const = 0; // is the current entry valid?
   virtual operator bool()       const = 0; //           ''

   virtual const T& operator *() const = 0; // return the current entry
   virtual T&       data()       const = 0; // return non-const entry
   virtual void update(const T&) const = 0; // update current entry

   virtual void     reset()            = 0; // return to first entry

   virtual bool     next(T&)           = 0; // return the next entry
   virtual bool     operator++()       = 0; // advance entry (prefix)
   virtual bool     operator++(int)    = 0; //    ''     ''  (postfix)
};

#endif
