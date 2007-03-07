// Name: Point3d.h
// Author: J. Michael Word
// Date Written: 10/28/96
// Purpose: Interface for Point3d and Matrix4d classes

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef POINT3D_H
#define POINT3D_H 1

#include <iostream>
#ifdef OLD_STD_HDRS
#include <math.h>
#else
#include <cmath>
using std::sqrt;
using std::sin;
using std::cos;
using std::acos;
#endif

// Note: all angles are in degrees. Conversion factors supplied below.

class Point3d;

typedef double Coord;
typedef Point3d Vector3d;

class Matrix4d { // homogenous coordinates transformation matrix
public:
   Matrix4d() { makeIdentity(); }               // identity matrix
   Matrix4d(Coord scale);                       // uniform scale matrix
   Matrix4d(Coord sx, Coord sy, Coord sz);      // scale matrix
   Matrix4d(const Point3d& point);              // translation matrix
   Matrix4d(const Vector3d& axis, double theta);// rotation matrix

   ~Matrix4d() { }

   Matrix4d operator*(const Matrix4d&) const; // matrix multiplication

   // to transform, we (post-)multiply vectors --> mat*vec
   Vector3d operator*(const Vector3d&) const;

   Matrix4d& makeIdentity();

   Coord _element[4][4]; // simple enough to be available to the public
};

// linear interpolation from lo (when a=0) to hi (when a=1)
inline Coord lerp(Coord lo, Coord hi, Coord a) {
   return lo + a*(hi-lo);
}

const Coord DEG2RAD=0.017453; // convert degrees to radians
const Coord RAD2DEG=57.29578; // convert radians to degrees

class Point3d { // point in R3
public:
   Point3d(): _x(0.0), _y(0.0), _z(0.0) { } // constructors
   Point3d(Coord x0, Coord y0, Coord z0): _x(x0),_y(y0),_z(z0){}
   Point3d(const Point3d& p): _x(p._x), _y(p._y), _z(p._z) { }
   ~Point3d() { }

   // create new points by adding or scaling
   Point3d operator-(const Point3d& p) const {
      return Point3d(_x - p._x, _y - p._y, _z - p._z);
   }
   Point3d operator+(const Point3d& p) const {
      return Point3d(_x + p._x, _y + p._y, _z + p._z);
   }
   Point3d operator*(Coord s) const {
      return Point3d(_x*s, _y*s, _z*s);
   }
   Point3d operator/(Coord s) const {
      return (s != 0.0) ?
               Point3d(_x/s, _y/s, _z/s) : Point3d(*this);
   }
   Point3d operator-() const {
      return Point3d(-_x, -_y, -_z);
   }

   Point3d& operator=(const Point3d& p) { // assignment
      _x = p._x; _y = p._y; _z = p._z;
      return *this;
   }
   Point3d& operator-=(const Point3d&); // update
   Point3d& operator+=(const Point3d&);
   Point3d& operator*=(Coord);
   Point3d& operator/=(Coord);

   Coord x() const { return _x; } // get
   Coord y() const { return _y; }
   Coord z() const { return _z; }
   const Point3d& loc() const { return *this; }

   Point3d& x(Coord x0) { _x = x0; return *this; } // set
   Point3d& y(Coord y0) { _y = y0; return *this; }
   Point3d& z(Coord z0) { _z = z0; return *this; }

   Coord lengthSquared() const {       // properties
      return (_x*_x) + (_y*_y) + (_z*_z);
   }
   Coord length() const {
      return sqrt(lengthSquared());
   }
   Vector3d& normalize(); // set length to unity
   Vector3d normal() const;

   Vector3d& scaleTo(Coord); // set length to a new value
   Vector3d scaled(Coord) const;

   // rotate (in degrees) our point around the vector(p1-p2)
   Point3d rotate(Coord, const Point3d&,const Point3d&) const;

private:
   Coord _x, _y, _z;
};

// dot product
inline Coord dot(const Point3d& a, const Point3d& b) {
      return (a.x()*b.x()) + (a.y()*b.y()) + (a.z()*b.z());
}
// cross product
Point3d cross(const Point3d&, const Point3d&);
// linear interpolation between two points or vectors
Point3d lerp(const Point3d&, const Point3d&, Coord);

// distance, etc. between two points
// apl - 2007/03/07 - inline for 10% speedup.
inline
Coord distanceSquared(const Point3d& a, const Point3d& b) {
   return Point3d( a - b ).lengthSquared();
}

Coord distance2(const Point3d&, const Point3d&);

// create a unit length vector in the pointing from b towards a
Vector3d makeVec(const Point3d& a, const Point3d& b);

// calculate angles and dihedrals in degrees
Coord angle(const Point3d&, const Point3d&, const Point3d&);
Coord dihedral(const Point3d&, const Point3d&,
               const Point3d&, const Point3d&);

// stream output
std::ostream& operator<<(std::ostream&, const Point3d&);
std::ostream& operator<<(std::ostream&, const Matrix4d&);
#endif
