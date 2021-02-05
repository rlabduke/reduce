// name: Rot3Fold.h
// author: J. Michael Word
// date written: 2/7/98
// purpose: Interface for Rot3Fold

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// Copyright (C) 2021 ReliaSolve LLC
// **************************************************************

#pragma once
#include "Rot.h"

/// @brief Insert and optimize a 3-fold Hydrogen set.
///
/// This includes methyl (CH3) groups and NH3 groups.  There is code to
/// fully optimize Methyl rotations, but it was found to be ineffective for
/// structures at the resolutions available up through 2020 and it seems like
/// much higher resolution structures would use something other than Reduce
/// to optimize these.  As a result, the Methyl optimization was removed from
/// Reduce and is not called.  The exception is Aromatic Methyl rotations, which
/// have 2 possible orientations.  They remain as a subclass of this method.
class Rot3Fold: public Rot {
public:
   Rot3Fold(const Point3d& a, const Point3d& b,
             const double ang, const PDBrec& heavyAtom);
   virtual ~Rot3Fold() {
   }

   virtual bool hasHires() const { return TRUE; }
   virtual void finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel, 
                         AtomPositions &xyz, DotSphManager& dotBucket);

   virtual int numOrientations(SearchStrategy ss=Mover::LOW_RES) const;
   virtual bool isDefaultO(int oi, SearchStrategy ss=Mover::LOW_RES) const;
   virtual std::string describeOrientation() const;

   virtual std::string formatComment(std::string prefix) const;

protected:
   virtual double orientationAngle(int oi, SearchStrategy ss=Mover::LOW_RES) const;
   virtual double orientationPenalty(float pmag) const;

   char        _grpName[20];

   Rot3Fold(const Rot3Fold& m); // copy and assign not implemented
   Rot3Fold& operator=(const Rot3Fold& m);
};
