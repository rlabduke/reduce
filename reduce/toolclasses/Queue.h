// Name: Queue.h
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: implement a object queue from a list
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#ifndef QUEUE_H
#define QUEUE_H 1

#include "List.h"
#include "ListIter.h"

template <class T>
class Queue : protected List<T> {
public:
   Queue() : List<T>() {}
   Queue(const Queue& q) : List<T>(q) {}

   bool empty() const { return List<T>::empty(); }
   int   size() const { return List<T>::size();  }

   void add(T obj) { insertBefore(obj, List<T>::END); }

   T remove() {
      T obj = List<T>::operator*();
      drop(List<T>::START);
      return obj;
   }
   const T& hd() const { return List<T>::operator*(); }

   friend ostream& operator << (ostream& s, const Queue<T>& l);
};

template <class T>
inline ostream& operator << (ostream& s, const Queue<T>& l) {
   return s << ((List<T>&)l);
}
#endif
