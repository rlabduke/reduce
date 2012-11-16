// name: DotSph.C
// author: J. Michael Word (port from dcr and mez fortran code)
// date written: 2/20/96 and converted to C++ 10/31/96
// purpose: Implementation for DotSphRep, DotSph and DotSphManager
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

#ifdef OLD_STD_HDRS
#include <math.h>
#else
#include <cmath>
using std::sin;
using std::cos;
using std::floor;
#endif

#include "DotSph.h"

// local utility functions
int estNumDots(float radius, float density);
int makeDots(float radius, Point3d points[], int maxpnts);

// make a dot sphere rep object including the dot array
DotSphRep::DotSphRep(float rad, float dens): _n(0), _p(0),
   _rad(rad), _dens(dens) {

   build(rad, dens);
}

// build a spherical shell of points
int DotSphRep::build(float radius, float density) {
   _rad  = radius;
   _dens = density;
   _n = 0;

   int m = estNumDots(_rad, _dens);
   _p = new Point3d[m + 1]; // extra point acts as our "nil" val
   if (_p) { _n = makeDots(_rad, &(_p[1]), m); }
   return _n;
}

// (slightly over)estimate of the number of dots
// makeDots() will use this dot count to determine the density
int estNumDots(float radius, float density) {
   return (int)floor(4.0 * PI * density * (radius * radius));
}

// generate 3D point coordinates at a specified radius
// with a density determined by maxpnts
int makeDots(float radius, Point3d points[], int maxpnts) {
	float offset = 0.2;
	double ang, cosang, sinang, phi, theta, xy0, x0, y0, z0;
	int i, j, k, odd, nequator, nvert, nhoriz;

	nequator = (int)floor(sqrt(maxpnts * PI));

	odd = 1;
	ang = 5.0 * PI / 360.0;
	cosang = cos(ang);
	sinang = sin(ang);

	i = 0;
	nvert = nequator / 2;
	for (j = 0; j <= nvert; j++) {
		phi = (PI * j) / nvert;
		z0 = cos(phi) * radius;
		xy0= sin(phi) * radius;

		nhoriz = (int)floor(nequator * sin(phi));
		if (nhoriz < 1) nhoriz = 1;
		for (k = 0; k < nhoriz; k++) {
			if(odd) {theta = (2.0 * PI * k + offset)/nhoriz; }
			else    {theta = (2.0 * PI * k         )/nhoriz; }
			x0 = cos(theta) * xy0;
			y0 = sin(theta) * xy0;

			if (i >= maxpnts) return i;
			points[i].x(x0);
			points[i].y(y0*cosang - z0*sinang);
			points[i].z(y0*sinang + z0*cosang);
			i++;
		}
		odd = !odd;
	}
	return i;
}

// copy constructor
DotSphManager::DotSphManager(const DotSphManager& m):
      _list(m._list), _dens(m._dens), _radFuzz(m._radFuzz),
      _densFuzz(m._densFuzz) {}
// assignment operator
DotSphManager& DotSphManager::operator=(const DotSphManager& m) {
   _list = m._list;
   _dens = m._dens;
   _radFuzz = m._radFuzz;
   _densFuzz = m._densFuzz;
   return *this;
}

DotSphManager::~DotSphManager()
{
   std::for_each(_list.begin(), _list.end(), DeleteObject());
}


// return a DotSph with a given radius & the default density
DotSph& DotSphManager::fetch(float rad) {
	DotSph *r;
	for(std::list<DotSph*>::const_iterator l = _list.begin(); l != _list.end(); ++l) {
		if (abs((*l)->radius() - rad)  < _radFuzz
			&& abs((*l)->density() - _dens) < _densFuzz) 
			return **l;
	}
	r = new DotSph();
	r->build(rad,_dens);
	_list.push_front(r);
	return *r;
}

// set the default dot density and return the previous value
float DotSphManager::density(float dens) {
   float od = _dens;
   _dens = dens;
   return od;
}

// set fuzziness used for comparing spheres and return prev val
float DotSphManager::radFuzz(float rf) {
   float orf = _radFuzz;
   _radFuzz = rf;
   return orf;
}
float DotSphManager::densFuzz(float df) {
   float odf = _densFuzz;
   _densFuzz = df;
   return odf;
}
