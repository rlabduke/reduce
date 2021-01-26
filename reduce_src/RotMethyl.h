// name: RotMethyl.h
// author: J. Michael Word
// date written: 2/7/98
// purpose: Interface for RotMethyl

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

class RotMethyl: public Rot {
public:
   RotMethyl(const Point3d& a, const Point3d& b,
             const double ang, const PDBrec& heavyAtom);
   virtual ~RotMethyl() {
   }

   virtual Mover::MemoType type() { return Mover::ROTATE_METHYL; }
   virtual bool hasHires() const { return TRUE; }
   virtual void finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel, 
                         AtomPositions &xyz, DotSphManager& dotBucket);
   virtual int makebumpers(std::multimap<LocBlk, std::shared_ptr<BumperPoint> >& bbins,
                           int n, float& maxVDWrad);
   virtual std::list<AtomDescr> getAtDescOfAllPos(float &maxVDWrad);

   virtual int numOrientations(SearchStrategy ss=Mover::LOW_RES) const;
   virtual bool isDefaultO(int oi, SearchStrategy ss=Mover::LOW_RES) const;
   virtual double orientationAngle(int oi, SearchStrategy ss=Mover::LOW_RES) const;
   virtual std::string describeOrientation() const;

   virtual std::string formatComment(std::string prefix) const;

protected:
   virtual bool setOrientation(int oi, float delta, AtomPositions &xyz,
	  SearchStrategy ss=Mover::LOW_RES);
   virtual double orientationPenalty(float pmag) const;

   char        _grpName[20];

   RotMethyl(const RotMethyl& m); // copy and assign not implemented
   RotMethyl& operator=(const RotMethyl& m);

};
