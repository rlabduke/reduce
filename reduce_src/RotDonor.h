// name: RotDonor.h
// author: J. Michael Word
// date written: 2/7/98
// purpose: Interface for RotDonor

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

class RotDonor: public Rot {
public:
   RotDonor(const Point3d& a, const Point3d& b,
            const double ang, const PDBrec& heavyAtom);
   virtual ~RotDonor() {};

   virtual Mover::MemoType type() { return Mover::ROTATE_DONOR; }
   virtual bool hasHires() const { return TRUE; }
   virtual void finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel, 
                         AtomPositions &xyz, DotSphManager& dotBucket);
   virtual int makebumpers(std::multimap<LocBlk, std::shared_ptr<BumperPoint> >& bbins,
                           int n, float& maxVDWrad);
   virtual std::list<AtomDescr> getAtDescOfAllPos(float &maxVDWrad);

   virtual int numOrientations(SearchStrategy ss=Mover::LOW_RES) const;
   virtual bool isDefaultO(int, SearchStrategy ss=Mover::LOW_RES) const { return FALSE; }
   virtual std::string describeOrientation() const;

   virtual std::string formatComment(std::string prefix) const;

protected:
   bool setOrientation(int oi, float delta, AtomPositions &xyz,
      SearchStrategy ss=Mover::LOW_RES);
   double orientationPenalty(float) const { return 0.0; }

   double bumpsThisAngle(AtomPositions &xyz, DotSphManager& dotBucket);

   int                   _nori;  // number of orientation angles
   std::vector< double >      _oang;  // array of angles

   RotDonor(const RotDonor& m); // copy and assign not implemented
   RotDonor& operator=(const RotDonor& m);
};
