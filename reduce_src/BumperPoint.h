// name: BumperPoint.h
// author: J. Michael Word
// date written: 2/7/98
// purpose: Interface for BumperPoint

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef BUMPERPOINT_H
#define BUMPERPOINT_H 1

#include "Point3d.h"
// #include "Hdl.h"
#include "utility.h"

class BumperPoint {
private:
   // internal representation of BumperPoint
   struct BumperPointRep {
      void update(const Point3d& p, int rn, int an, float rad) {
	 _loc=p; _rnum=rn; _anum=an;
	 _radii=rad; _ok=TRUE;
      }
      Point3d _loc;
      int     _rnum;
      int     _anum;
      float   _radii;
      bool    _ok;
   };

	//private and unimplemented
	BumperPoint(const BumperPoint& b); //: _rep(b._rep) {}
	BumperPoint& operator=(const BumperPoint& b);
	//{
   //	_rep = b._rep;
   //	return *this;
   //}

public:
   BumperPoint(const Point3d& p, int rn, int an, float rad) {
	   _rep = new BumperPointRep();
	   _rep->update(p,rn,an,rad); 
   }
  ~BumperPoint() { delete _rep;}
   void invalidate() { _rep->_ok = FALSE; };

   const Point3d& loc() const { return _rep->_loc; }

   int     rnum() const { return _rep->_rnum; };
   int     anum() const { return _rep->_anum; };
   float vdwRad() const { return _rep->_radii; };
   bool   valid() const { return _rep->_ok; };
private:
   BumperPointRep* _rep;
};
#endif
