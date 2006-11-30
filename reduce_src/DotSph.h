// name: DotSph.C
// author: J. Michael Word (port from dcr and mez fortran code)
// date written: 2/20/96 and converted to C++ 10/31/96
// purpose: Interface for DotSphRep, DotSph and DotSphManager
//          which create and manage spherical arrays of points

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#if defined(_MSC_VER)
#pragma warning(disable:4305) 
#endif

#ifndef DOTSPH_H
#define DOTSPH_H 1

#include "utility.h"
#include "Point3d.h"
#include <list>
#include <algorithm>

// Because the list of points in a dot sphere can be large
// the structure has been broken into two classes. The "rep"
// class is shared rather than ever being actually copied
// and the DotSph class is what the user actually works with.

// representation class for spherical shells of dots 
class DotSphRep {
public:
   DotSphRep(): _n(0), _p(0), _rad(0), _dens(0) {}
   DotSphRep(float rad, float dens);
   ~DotSphRep() { if (_p) delete [] _p; }

   DotSphRep(const DotSphRep&); // unimplemented to prevent copy
   DotSphRep& operator=(const DotSphRep&); // unimplemented to prevent copy

   bool operator==(const DotSphRep& d) { //comparison operators
      return _rad == d._rad && _dens == d._dens;
   }
   bool operator!=(const DotSphRep& d) {
      return _rad != d._rad || _dens != d._dens;
   }
   int build(float rad, float dens);
   int count() const { return _n; }
   // note: internally we use the point 0 as a "nil" value
   const Point3d* points() const { return _p; }
   float radius()  const { return _rad; }
   float density() const { return _dens; }
private:
   int      _n;          // number of three-dimensional points
   Point3d *_p;          // array of points

   float _rad;            // sphere radius
   float _dens;           // dot density
};

// visible interface for dot spheres
class DotSph {
public:
	DotSph() { _rep = NULL; }
	DotSph(float rad, float dens) { build(rad, dens); }
   ~DotSph() {
	   delete _rep; //delete 0 is fine
   }

private:
	//apl private and unimplemented since deleting rep without looking at a 
	//apl counter (as is done in class PDB with the int* PDB::_i ) makes
	//apl dangling references likely
   DotSph(const DotSph& s); //{ _rep = s._rep; }
   DotSph& operator=(const DotSph& s);
	//{
   //   _rep = s._rep; return *this;
   //}
public:

   // strict comparisons
   bool operator==(const DotSph& s) { return *(_rep) == *(s._rep); }
   bool operator!=(const DotSph& s) { return *(_rep) != *(s._rep); }

   int build(float rad, float dens) { // forget old dots and rebuild
      _rep = new DotSphRep();     // don't affect any other copies!
      return _rep->build(rad, dens);
   }
   int count() const { return _rep->count(); } // how many dots?

    // return the ith dot in the range [0..count-1] or the nil dot
    // note: we have to account for rep's use of the 0th dot
   const Point3d& operator[](int i) const {
      return (_rep->points())[(i < 0 || i >= count()) ? 0 : i+1];
   }

   // return the sphere parameters
   float radius()  const { return _rep->radius(); }
   float density() const { return _rep->density(); }
private:
	DotSphRep *_rep;
};

// manage a cache of dot spheres, creating or reusing them as required
class DotSphManager {
public:
   DotSphManager(): _dens(16.0), _radFuzz(0.001), _densFuzz(0.1) {}
   DotSphManager(float d): _dens(d), _radFuzz(0.001), _densFuzz(0.1) {}
   ~DotSphManager();
	   
   DotSphManager(const DotSphManager& m);
   DotSphManager& operator=(const DotSphManager& m);

   DotSph& fetch(float rad);
   int count() const { return _list.size(); } // how many spheres?

   float density(float); // set the default dot density

   float radFuzz(float); // set the comparison fuzz factors
   float densFuzz(float);

private:
	std::list<DotSph*> _list; // list of spheres already built
   float _dens;       // default dot density per square angstrom

   float _radFuzz;    // fuzziness when comparing spheres
   float _densFuzz;
};

#endif
