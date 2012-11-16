// Name: ListIter.C
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: List iterator implementation
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#include "ListIter.h"

template <class T>
bool ListIter<T>::next(T& e) {   // return current entry and advance
   bool rc = (_node != NULL);
   if (rc) {
      e = _list.linkData(_node);
      _node = _list.linkNext(_node);
   }
   return rc;
}

template <class T>
bool ListIter<T>::prev(T& e) {   // return current entry and retreat
   bool rc = (_node != NULL);
   if (rc) {
      e = _list.linkData(_node);
      _node = _list.linkPrev(_node);
   }
   return rc;
}

template <class T>
bool ListIter<T>::operator++() {    // advance to next entry (prefix)
   assert(_node);
   _node = _list.linkNext(_node);
   return _node != NULL;
}

template <class T>
bool ListIter<T>::operator++(int) { // advance to next entry (postfix)
   assert(_node);
   _node = _list.linkNext(_node);
   return _node != NULL;
}

template <class T>
bool ListIter<T>::operator--() {    // retreat to prev entry (prefix)
   assert(_node);
   _node = _list.linkPrev(_node);
   return _node != NULL;
}

template <class T>
bool ListIter<T>::operator--(int) { // retreat to prev entry (postfix)
   assert(_node);
   _node = _list.linkPrev(_node);
   return _node != NULL;
}

template <class T>
bool NonConstListIter<T>::drop() {
   bool rc = (_node != NULL);
   if (rc) {
      DblLnkLstNode<T>* prev = _NClist.linkPrev(_node);
      _NClist.drop(_node);
      _node = prev;
   }
   return rc;
}

template <class T>
bool NonConstListIter<T>::insertBefore(const T& e) {
   bool rc = (_node != NULL);
   if (rc) {
      _NClist.insertBefore(e, _node);
   }
   return rc;
}

template <class T>
bool NonConstListIter<T>::insertAfter(const T& e) {
   bool rc = (_node != NULL);
   if (rc) {
      _NClist.insertAfter(e, _node);
   }
   return rc;
}

