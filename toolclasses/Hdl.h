// Name: Hdl.h
// Author: J. Michael Word
// Date Written: 10/25/96
// Purpose: Interface for Hdl class
//          which are used to manage a surrogate of an object.
//          From: Ruminations on C++, Koenig and Moo, Addison Wesley, 1996
#ifndef HDL_H
#define HDL_H 1
#include "UseCount.h"

#ifndef BOOLPREDEFINED
typedef int bool;
#endif
#ifdef BRACKETOPERPARMS
#define PARM_TYPE_T <T>
#else
#define PARM_TYPE_T
#endif

template <class T> class Hdl {
public:
   Hdl(): _p(new T) { }
   Hdl(const T& t): _p(new T(t)) { }
   Hdl(T* t): _p(t) { }
   Hdl(const Hdl& h): _u(h._u), _p(h._p) { }
   ~Hdl() { if (_u.only()) delete _p; }

   Hdl& operator=(const Hdl&);

   friend bool operator== PARM_TYPE_T(const Hdl&, const Hdl&);
   friend bool operator!= PARM_TYPE_T(const Hdl&, const Hdl&);
   friend bool operator<  PARM_TYPE_T(const Hdl&, const Hdl&);
   friend bool operator>  PARM_TYPE_T(const Hdl&, const Hdl&);
   friend bool operator<= PARM_TYPE_T(const Hdl&, const Hdl&);
   friend bool operator>= PARM_TYPE_T(const Hdl&, const Hdl&);

   const T& operator*() const { return *_p; }
   T& data() { return *_p; }

   void makeonly() { if (_u.makeonly()) { _p = new T(*_p); } }
private:
   T* _p;
   UseCount _u;
};

// I have made makeonly() and data() public so that derrived
// and containing classes may write functions to access and
// set T variables.
//
//   PointHandle& PointHandle::x(int xx) {
//      _u.makeonly(); <--- if semantics require a copy on write
//      data()._x = xx;
//      return *this;
//   }

template <class T>
inline bool operator==(const Hdl<T>& op1, const Hdl<T>& op2) {
   return (*(op1._p) == *(op2._p));
}
template <class T>
inline bool operator!=(const Hdl<T>& op1, const Hdl<T>& op2) {
   return !(*(op1._p) == *(op2._p));
}
template <class T>
inline bool operator<(const Hdl<T>& op1, const Hdl<T>& op2) {
   return (*(op1._p) < *(op2._p));
}
template <class T>
inline bool operator>(const Hdl<T>& op1, const Hdl<T>& op2) {
   return (*(op1._p) > *(op2._p));
}
template <class T>
inline bool operator<=(const Hdl<T>& op1, const Hdl<T>& op2) {
   return (*(op1._p) <= *(op2._p));
}
template <class T>
inline bool operator>=(const Hdl<T>& op1, const Hdl<T>& op2) {
   return (*(op1._p) >= *(op2._p));
}

#ifdef INCTEMPLATEDEFNS
#include "Hdl.cpp"
#endif

#endif
