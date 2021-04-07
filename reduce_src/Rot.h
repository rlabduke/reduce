// name: Rot.h
// author: J. Michael Word
// date written: 2/7/98
// purpose: Base class for all rotation Movers

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

#include "PDBrec.h"
#include "DotSph.h"
#include "BumperPoint.h"
#include "Mover.h"
#include <vector>
#include "neighbors.h"
#include "utility.h"

class Rot: public Mover {
public:
  Rot(const Point3d& a, const Point3d& b,
            const double ang, const PDBrec& heavyAtom);
  virtual ~Rot() {
  }

  virtual bool isComplete() const { return TRUE; }
  virtual bool canFlip() const { return FALSE; }
  virtual bool canRotate() const { return TRUE; }
  virtual int flipState() const { return 0; }
  virtual bool markFlipAtoms() { return FALSE; }
  virtual const PDBrec& exampleAtom() const { return heavyAtom(); }
  virtual void limitOrientations(bool, SearchStrategy ss=Mover::LOW_RES) { /* Do nothing */};
  virtual double determineScore(AtomPositions &xyz,
    DotSphManager& dotBucket, int nBondCutoff,
    float probeRadius, float pmag, double& penalty,
    float &bumpScore, float &hbScore, bool& hasBadBump);

  virtual int makebumpers(std::multimap<LocBlk, std::shared_ptr<BumperPoint> >& bbins,
                          int n, float& maxVDWrad);
  virtual std::list<AtomDescr> getAtDescOfAllPos(float &maxVDWrad);

  const PDBrec& heavyAtom() const { return *_heavyAtom; }

  virtual bool insertHAtom(const PDBrec& ha) {
      std::shared_ptr<PDBrec> temp = std::make_shared<PDBrec>();
    *temp = ha;
    _rot.push_front(temp);
    return true;
  }

  virtual bool setOrientation(int oi, AtomPositions &xyz, SearchStrategy ss=Mover::LOW_RES) {
     return setOrientation(oi, 0.0, xyz, ss);
  }

  double angle() const { return _angle; }

  virtual void setHydAngle(double newAng, AtomPositions &xyz);
  virtual void dropBondedFromBumpingListForPDBrec( std::list< std::shared_ptr<PDBrec> > & bumping, std::shared_ptr<PDBrec> atom, int nBondCutoff  ) const;

protected:
  virtual double orientationAngle(int oi, SearchStrategy ss=Mover::LOW_RES) const = 0;
  virtual bool setOrientation(int oi, float delta, AtomPositions &xyz,
	  SearchStrategy ss=Mover::LOW_RES);
  void angle(double val) { _angle = clampAngle(val); }
  double scoreThisAngle(AtomPositions &xyz,
     DotSphManager& dotBucket, int nBondCutoff,
     float probeRadius, float &bumpScore,
     float &hbScore, bool &hasBadBump);
  int findAtom(std::shared_ptr<PDBrec> atom ) const;
  virtual void setHydAngle(PDBrec& theAtom, double oldAng, double newAng,
                        AtomPositions &xyz);

  Point3d     _p1, _p2;   // rotation axis is through these points
  std::shared_ptr<PDBrec>      _heavyAtom = std::make_shared<PDBrec>(); // hydrogen attachment point
  std::list< std::shared_ptr<PDBrec> > _rot;       // rotating hydrogen atoms
  double      _angle;

  std::vector< std::shared_ptr<std::list< std::shared_ptr<PDBrec> > > > _bnded; // pre-calculated bonding list

  Rot(const Rot& m); // copy and assign not implemented
  Rot& operator=(const Rot& m);

  // These change the way the orientation methods work.
  // Derived classes have different values from the base class.
  double START_ANGLE = 180;
  double ROUGH_STEP = 30;
  double FINE_STEP = 1;

  // These change the way the makebumbers() method behaves and are
  // modified in the constructors of derived classes
	double dtheta;   ///< fineness of rotation angle scan
	double scanAngle;
};
