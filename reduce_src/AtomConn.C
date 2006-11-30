// name: AtomConn.C
// author: J. Michael Word
// date written: 7/15/97
// purpose: Connected atoms and proton placement plans

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
#pragma warning(disable:4786) 
#endif

#include "AtomConn.h"

AtomConn& AtomConn::operator=(const AtomConn& a) {
   if (this != &a) {
      _name     = a._name;
      _order    = a._order;
      _neighbor = a._neighbor;
	  _emptyString = "";
   }
   return *this;
}

// return the "i"th connection [0..num_conn-1] for this atom
const std::string& AtomConn::conn(int i) const {
	if (i >= num_conn() || i < 0) { return _emptyString; }
   return _neighbor[i];
}
   
void AtomConn::limitConnections(int n) {
//   int n2cut = (n >= 0) ? (num_conn() - n) : 0;
	if (n >= 0)
		_neighbor.erase(_neighbor.begin() + n, _neighbor.end());
}

atomPlacementPlan::atomPlacementPlan(const atomPlacementPlan& a)
      : AtomConn(a), _type(a._type), _elem(a._elem), _dist(a._dist),
                     _ang1(a._ang1), _ang2(a._ang2), _flags(a._flags) {}

atomPlacementPlan& atomPlacementPlan::operator=(const atomPlacementPlan& a) {
   if (this != &a) {
      AtomConn::operator=(a);
      _type  = a._type;
      _elem  = a._elem;
      _dist  = a._dist;
      _ang1  = a._ang1;
      _ang2  = a._ang2;
      _flags = a._flags;
   }
   return *this;
}

bool atomPlacementPlan::placeH(const std::vector<Point3d>& loc, Point3d& hpos) const {
   bool rc = FALSE;

   switch(_type) {
   case 1:
      if (loc.size() == 4) {
	 hpos = type1position(loc[0], loc[1], loc[2], loc[3], _dist);
	 rc = TRUE;
      }
      break;
   case 2:
      if (loc.size() == 3) { // in this case ang2 is really a fudge factor
	 hpos = type2position(loc[0], loc[1], loc[2], _dist, _ang1, _ang2);
	 rc = TRUE;
      }
      break;
   case 3:
      if (loc.size() == 3) {
	 hpos = type3position(loc[0], loc[1], loc[2], _dist, _ang1, _ang2);
	 rc = TRUE;
      }
      break;
   case 4:
      if (loc.size() == 3) { // in this case ang2 is really a fudge factor
	 hpos = type4position(loc[0], loc[1], loc[2], _dist, _ang2);
	 rc = TRUE;
      }
      break;
   case 5:
      if (loc.size() == 3) { // in this case ang2 is really a fraction
	 hpos = type5position(loc[0], loc[1], loc[2], _dist, _ang2);
	 rc = TRUE;
      }
      break;
   case 6:
      if (loc.size() == 2) {
	 hpos = type6position(loc[0], loc[1], _dist);
	 rc = TRUE;
      }
      break;
   }
	
	//if (!rc) std::cerr << "Failed to place hydrogen: " << _type << " " << loc.size() << std::endl;
	
   return rc;
}

Point3d atomPlacementPlan::calcLoc(int locType,
                      const Point3d& a1, const Point3d& a2,
                      const Point3d& a3, const Point3d& a4,
                      float len, float ang, float xtra) {
   Point3d hpos(-999.9,-999.9,-999.9);

   switch(locType) {
   case 1: hpos = type1position(a1, a2, a3, a4, len);            break;
   case 2: hpos = type2position(a1, a2, a3,     len, ang, xtra); break;
   case 3: hpos = type3position(a1, a2, a3,     len, ang, xtra); break;
   case 4: hpos = type4position(a1, a2, a3,     len,      xtra); break;
   case 5: hpos = type5position(a1, a2, a3,     len,      xtra); break;
   case 6: hpos = type6position(a1, a2,         len);            break;
   }
   return hpos;
}

// 1: HXR3 - requires just 4 atom centers
Point3d atomPlacementPlan::type1position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos, const Point3d& a4pos,
                      float bondlen) {

   const Vector3d v12 = (a2pos - a1pos).normalize();
   const Vector3d v13 = (a3pos - a1pos).normalize();
   const Vector3d v14 = (a4pos - a1pos).normalize();

   const Vector3d bondvec = (v12 + v13 + v14).scaleTo(-bondlen);

   return a1pos + bondvec;
}

// 2: H2XR2- three atoms and an angle
Point3d atomPlacementPlan::type2position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos,
                      float bondlen, float angle, float fudge) {

   const Vector3d v12 = (a2pos - a1pos).normalize();
   const Vector3d v13 = (a3pos - a1pos).normalize();

   const Point3d between = a1pos + lerp(v12, v13, 0.5 + fudge);

   return atomPlacementPlan::type3position(
                      a1pos, between, a2pos, bondlen, angle, 90.0);
}

// 3: H3XR - three atoms an angle and dihedral
Point3d atomPlacementPlan::type3position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos,
                      float bondlen, float theta, float phi) {

   const Vector3d  v21 = (a1pos - a2pos).normalize();
   const Vector3d  v23 = (a3pos - a2pos).normalize();
   const Vector3d norm = cross(v21, v23).scaleTo(bondlen);

   const Point3d pos4 = (a1pos + norm).rotate(phi-90.0, a2pos, a1pos);

   const Vector3d v14 = (pos4  - a1pos).normalize();
   const Vector3d v12 = (a2pos - a1pos).normalize();
   const Point3d pos5 = cross(v14, v12) + a1pos;
 
   return pos4.rotate(90.0-theta, a1pos, pos5);
}

// 4: HXR2 - three atoms and a fudge factor
Point3d atomPlacementPlan::type4position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos,
                      float bondlen, float fudge) {

   const Vector3d v12 = (a2pos - a1pos).normalize();
   const Vector3d v13 = (a3pos - a1pos).normalize();

   const Vector3d hvec = lerp(v12, v13, 0.5 + fudge).scaleTo(-bondlen);

   return a1pos + hvec;
}

// 5: HXR2 - three atoms and a fraction
Point3d atomPlacementPlan::type5position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos,
                      float bondlen, float fract) {

   const Vector3d v12 = (a2pos - a1pos).scaleTo(bondlen);
   const Vector3d v13 = (a3pos - a1pos).scaleTo(bondlen);
   const Point3d pos4 = cross(v12, v13) + a1pos;

   const Coord cncaAngle = angle(a2pos, a1pos, a3pos);
   const Coord hncaAngle = fract*(360.0 - cncaAngle);

   return (a1pos + v12).rotate(hncaAngle, pos4, a1pos);
}

// 6: HXY  - (linear) just two atoms
Point3d atomPlacementPlan::type6position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      float bondlen) {

   const Vector3d hvec = (a1pos - a2pos).scaleTo(bondlen);

   return a1pos + hvec;
}
