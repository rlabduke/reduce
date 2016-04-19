// name: LocBlk.h
// author: J. Michael Word
// date written: 8/20/97
// purpose: Interface for LocBlk

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef LOCBLK_H
#define LOCBLK_H 1

#ifdef OLD_STD_HDRS
#include <math.h>
#else
#include <cmath>
using std::fabs;
using std::floor;
#endif

#include "Point3d.h"
#include "utility.h"

const Coord LocBlkTolerance = 0.00001;
const Coord LocBlkGrain     = 0.3125;    // i.e. 3.2x3.2 Angstrom blocks

class LocBlk {
public:
   LocBlk() {}
   LocBlk(const LocBlk& b): _rep(b._rep) {}
   LocBlk(Coord x0, Coord y0, Coord z0) { chunk(_rep, Point3d(x0, y0, z0)); }
   LocBlk(const Point3d& p) { chunk(_rep, p); }
  ~LocBlk() {}

   LocBlk& operator=(const LocBlk& b) { _rep=b._rep; return *this; }
   LocBlk& operator=(const Point3d& p) {
      chunk(_rep, p);
      return *this;
   }

   Coord x() const { return _rep.x(); }
   Coord y() const { return _rep.y(); }
   Coord z() const { return _rep.z(); }
   const Point3d& loc() const { return _rep; }

   bool operator == (const LocBlk& b) const {
      return ::fabs(_rep.x() - b._rep.x()) < LocBlkTolerance
          && ::fabs(_rep.y() - b._rep.y()) < LocBlkTolerance
          && ::fabs(_rep.z() - b._rep.z()) < LocBlkTolerance;
   }
   bool operator < (const LocBlk& b) const {
      return _rep.x() < b._rep.x()
//          || (_rep.y() < b._rep.y())
//          || (_rep.z() < b._rep.z());
          || (_rep.x() == b._rep.x() && _rep.y() < b._rep.y())
          || (_rep.x() == b._rep.x() && _rep.y() == b._rep.y() && _rep.z() < b._rep.z());
   }

   static Coord chunksz(Coord c) { return c/LocBlkGrain; }
   static Coord scale(Coord c) { return c*LocBlkGrain; }

   static void chunk(Point3d& r, const Point3d& p) {
      r.x(::floor(p.x()*LocBlkGrain));
      r.y(::floor(p.y()*LocBlkGrain));
      r.z(::floor(p.z()*LocBlkGrain));
   }
private:

   Point3d _rep;
};

inline unsigned long hash(const LocBlk& key, unsigned long M) {
   union { long sl; unsigned long ul; };
   unsigned long h = 0;
   sl = long(key.x()*250000.0); h =      (ul)  % M;
   sl = long(key.y()*500.0);    h = (h + (ul)) % M;
   sl = long(key.z());          h = (h + (ul)) % M;
   return h;
}

#endif
