// Name: Hdl.h
// Author: J. Michael Word
// Date Written: 10/25/96
// Purpose: Implementation for Hdl class
//          which is used to manage a surrogate of an object.
//          From: Ruminations on C++, Koenig and Moo, Addison Wesley, 1996
#include "Hdl.h"

template <class T>
Hdl<T>& Hdl<T>::operator=(const Hdl& h) {
   if (_u.reattach(h._u)) { delete _p; }
   _p = h._p;
   return *this;
}
