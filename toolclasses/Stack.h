// Name: Stack.h
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: implement a object stack from a list
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#ifndef STACK_H
#define STACK_H 1

#include "List.h"
#include "ListIter.h"

template <class T>
class Stack : protected List<T> {
public:
   Stack() : List<T>() {}
   Stack(const Stack& s) : List<T>(s) {}

   bool empty() const { return List<T>::empty(); }
   int   size() const { return List<T>::size();  }

   void push(T obj) { insertBefore(obj, List<T>::START); }

   T pop() {
      T obj = List<T>::operator*();
      drop(List<T>::START);
      return obj;
   }
   const T& top() const { return List<T>::operator*(); }

   friend ostream& operator << (ostream& s, const Stack<T>& l);
};

template <class T>
inline ostream& operator << (ostream& s, const Stack<T>& l) {
   return s << ((List<T>&)l);
}
#endif
