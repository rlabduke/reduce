// Name: Point3d.C
// Author: J. Michael Word
// Date Written: 10/28/96
// Purpose: Implementation for Point3d and Matrix4d classes

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#include "Point3d.h"

// scale a vector to unit length
Vector3d& Point3d::normalize() {
   return Point3d::operator/=(length());
}

// generate a vector from a prev one scaled to unit length
Vector3d Point3d::normal() const {
   return Point3d::operator/(length());
}

// scale a vector to a given length
Vector3d& Point3d::scaleTo(Coord newlen) {
	return normalize() *= newlen;
}
// generate a vector from a prev one scaled to a given length
Vector3d Point3d::scaled(Coord newlen) const {
	return normal() * newlen;
}

Point3d& Point3d::operator-=(const Point3d& p) {
   _x -= p._x;
   _y -= p._y;
   _z -= p._z;
   return *this;
}
Point3d& Point3d::operator+=(const Point3d& p) {
   _x += p._x;
   _y += p._y;
   _z += p._z;
   return *this;
}
Point3d& Point3d::operator*=(Coord s) {
   _x *= s;
   _y *= s;
   _z *= s;
   return *this;
}
Point3d& Point3d::operator/=(Coord s) {
   if (s != 0.0) {
      _x /= s;
      _y /= s;
      _z /= s;
   }
   return *this;
}

Point3d cross(const Point3d& a, const Point3d& b) {
   return Point3d( (a.y()*b.z()) - (a.z()*b.y()),
	               (a.z()*b.x()) - (a.x()*b.z()),
	               (a.x()*b.y()) - (a.y()*b.x()) );
}

Point3d lerp(const Point3d& lo, const Point3d& hi, Coord alpha) {

   return Point3d( lerp(lo.x(), hi.x(), alpha),
                   lerp(lo.y(), hi.y(), alpha),
                   lerp(lo.z(), hi.z(), alpha) );
}


Coord distance2(const Point3d& a, const Point3d& b) {
	return (sqrt(distanceSquared(a, b)));
}

Vector3d makeVec(const Point3d& a, const Point3d& b) {
   return Point3d( a - b ).normalize();
}

// ---------- matrix related

Matrix4d::Matrix4d(Coord s) { // uniform scale matrix
   makeIdentity();
   _element[0][0] = s;
   _element[1][1] = s;
   _element[2][2] = s;
}

Matrix4d::Matrix4d(Coord sx, Coord sy, Coord sz) { // scale matrix
   makeIdentity();
   _element[0][0] = sx;
   _element[1][1] = sy;
   _element[2][2] = sz;
}

Matrix4d::Matrix4d(const Point3d& p) { // translation matrix
   makeIdentity();
   _element[0][3] = p.x();
   _element[1][3] = p.y();
   _element[2][3] = p.z();
}

// build a theta (degrees) angle rotation matrix with vector as an axis

Matrix4d::Matrix4d(const Vector3d& v, double theta) {
	const Coord c = cos(theta*DEG2RAD);
	const Coord s = sin(theta*DEG2RAD);
	const Coord t = 1.0 - c;
	const Coord x = v.x();
	const Coord y = v.y();
	const Coord z = v.z();

	makeIdentity();

	_element[0][0] = t * x * x    +        c;
	_element[0][1] = t * x * y    +    z * s;
	_element[0][2] = t * x * z    -    y * s;

	_element[1][0] = t * y * x    -    z * s;
	_element[1][1] = t * y * y    +        c;
	_element[1][2] = t * y * z    +    x * s;

	_element[2][0] = t * z * x    +    y * s;
	_element[2][1] = t * z * y    -    x * s;
	_element[2][2] = t * z * z    +        c;
}

Matrix4d& Matrix4d::makeIdentity() { // build an identity matrix
   _element[0][0] = 1.0;
   _element[0][1] = 0.0;
   _element[0][2] = 0.0;
   _element[0][3] = 0.0;

   _element[1][0] = 0.0;
   _element[1][1] = 1.0;
   _element[1][2] = 0.0;
   _element[1][3] = 0.0;

   _element[2][0] = 0.0;
   _element[2][1] = 0.0;
   _element[2][2] = 1.0;
   _element[2][3] = 0.0;

   _element[3][0] = 0.0;
   _element[3][1] = 0.0;
   _element[3][2] = 0.0;
   _element[3][3] = 1.0;

   return *this;
}

Vector3d Matrix4d::operator*(const Vector3d& v) const {
   Coord x, y, z, w;
   x = (v.x() * _element[0][0]) +
       (v.y() * _element[1][0]) +
       (v.z() * _element[2][0]) +
                _element[3][0]; 
   y = (v.x() * _element[0][1]) +
       (v.y() * _element[1][1]) +
       (v.z() * _element[2][1]) +
                _element[3][1]; 
   z = (v.x() * _element[0][2]) +
       (v.y() * _element[1][2]) +
       (v.z() * _element[2][2]) +
                _element[3][2]; 
   w = (v.x() * _element[0][3]) +
       (v.y() * _element[1][3]) +
       (v.z() * _element[2][3]) +
                _element[3][3]; 

   return Vector3d(x, y, z) /= w;
}

Matrix4d Matrix4d::operator*(const Matrix4d& m2) const {
   Matrix4d c;
   int i, j, k;
   for (i=0; i < 4; i++) {
      for (j=0; j < 4; j++) {
         c._element[i][j] = 0.0;
         for (k=0; k < 4; k++) {
            c._element[i][j] +=
                   _element[i][k] * m2._element[k][j];
         }
      }
   }
   return c;
}

// ---------- angles and dihedrals

// rotate point theta degrees around the a->b axis to yield a new point

Point3d Point3d::rotate(Coord theta,
                        const Point3d& a, const Point3d& b) const {
	Matrix4d rotmat(makeVec(b, a), theta);

	return (rotmat * Point3d::operator-(b)) + b;
}

// calculate the angle (degrees) between 3 points

Coord angle(const Point3d& p1, const Point3d& p2, const Point3d& p3) {
	Vector3d a = p1 - p2;
	Vector3d b = p3 - p2;

	Coord amag   = a.length();
	Coord bmag   = b.length();
	Coord theta  = 0.0;

	if (amag*bmag >= 0.0001) { theta = acos(dot(a, b)/(amag*bmag)); }

	return theta*RAD2DEG;
}

// calculate the dihedral angle (degrees) given 4 points

Coord dihedral(const Point3d& p1, const Point3d& p2,
               const Point3d& p3, const Point3d& p4) {
	Vector3d b; // used to set handedness

	Vector3d d = cross((    p1 - p2), (p3 - p2));
	Vector3d e = cross((b = p2 - p3), (p4 - p3));

	Coord dmag  = d.length();
	Coord emag  = e.length();
	Coord theta = 0.0;

	if (dmag*emag >= 0.0001) { theta = acos(dot(d, e)/(dmag*emag)); }

	Vector3d f = cross(d, b); // this part sets the correct handedness

	Coord fmag = f.length();
	Coord phi  = 0.0;

	if (fmag*emag >= 0.0001) { phi = acos(dot(f, e)/(fmag*emag)); }

	if (phi*RAD2DEG > 90.0) { theta = - theta; };

	return theta*RAD2DEG;
}

// ----- stream output

std::ostream& operator<<(std::ostream& os, const Point3d& p) {
   os << "{" << p.x() << ", " << p.y() << ", " << p.z() << "}";
   return os;
}

std::ostream& operator<<(std::ostream& os, const Matrix4d& m) {
   os << "{{" << m._element[0][0] << ", " << m._element[0][1] << ", "
              << m._element[0][2] << ", " << m._element[0][3] << "},";
   os << " {" << m._element[1][0] << ", " << m._element[1][1] << ", "
              << m._element[1][2] << ", " << m._element[1][3] << "},";
   os << " {" << m._element[2][0] << ", " << m._element[2][1] << ", "
              << m._element[2][2] << ", " << m._element[2][3] << "},";
   os << " {" << m._element[3][0] << ", " << m._element[3][1] << ", "
              << m._element[3][2] << ", " << m._element[3][3] << "}}";
   return os;
}
