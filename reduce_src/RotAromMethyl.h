// name: RotAromMethyl.h
// author: J. Michael Word, modified by Aram Han
// date written: 2/7/98, modified 8/13/12
// purpose: Interface for RotAromMethyl

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef ROTAROMMETHYL_H
#define ROTAROMMETHYL_H 1

#include "RotMethyl.h"

class RotAromMethyl: public RotMethyl {
public:
   RotAromMethyl(const Point3d& a, const Point3d& b,
             const double ang, const PDBrec& heavyAtom);
   virtual ~RotAromMethyl() { }

   virtual bool hasHires() const { return FALSE; }
   virtual int numOrientations(SearchStrategy ss=Mover::LOW_RES) const;

   virtual bool setOrientation(int oi, float delta, AtomPositions &xyz,
	                SearchStrategy ss=Mover::LOW_RES);
   
private:
   RotAromMethyl(const RotAromMethyl& m); // copy and assign not implemented
   RotAromMethyl& operator=(const RotAromMethyl& m);
};

#endif
