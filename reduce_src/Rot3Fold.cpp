// name: Rot3Fold.C
// author: J. Michael Word
// date written: 2/7/98
// purpose: Implementation for Rot3Fold

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

#ifdef OLD_STD_HDRS
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#else
#include <cstdio>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
using std::cerr;
using std::endl;
using std::strcpy;
using std::exit;
#endif

#include "Rot3Fold.h"
#include "AtomPositions.h"

Rot3Fold::Rot3Fold(const Point3d& a, const Point3d& b,
                     const double ang, const PDBrec& heavyAtom)
   : Rot(a, b, ang, heavyAtom)
{
  strcpy(_grpName, ((heavyAtom.elem().atno() == 7) ?
		"NH3+   " : "methyl "));

  // Specify values for this subclass
  START_ANGLE = 180;
  ROUGH_STEP = 30;
  FINE_STEP = 1;
  dtheta = 10;
  scanAngle = 60;
}

void Rot3Fold::finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel,
						 AtomPositions &xyz, DotSphManager& dotBucket) {

	if (isComplete()) {

		// pre-build lists of bonded atoms

		const double approxNbondDistLimit = 3.0 + 0.5*nBondCutoff; // just a rule of thumb

		_rot.push_front(_heavyAtom);
		for(std::list< std::shared_ptr<PDBrec> >::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
			const  std::shared_ptr<PDBrec> thisAtom = *alst;
			if (thisAtom->valid()) {
				std::shared_ptr<std::list< std::shared_ptr<PDBrec> > > temp = std::make_shared<std::list< std::shared_ptr<PDBrec> > >();
			  std::list< std::shared_ptr<PDBrec> > neighborList = xyz.neighbors(thisAtom->loc(), thisAtom->covRad(),
				approxNbondDistLimit);
			  bondedList(*thisAtom, neighborList, nBondCutoff, _rot, temp.get());
			  _bnded.push_back(temp);
            }
		}
		_rot.pop_front();
	}
}

int Rot3Fold::numOrientations(SearchStrategy ss) const {
   return (ss==Mover::LOW_RES) ?
      int(120.0/ROUGH_STEP + 0.5) :
      int(2.0*(double(ROUGH_STEP)/double(FINE_STEP)) + 0.5);
}

bool Rot3Fold::isDefaultO(int oi, SearchStrategy ss) const {
   if (ss!=Mover::LOW_RES) { oi = orientation(Mover::LOW_RES); }
   return ((oi == 0) &&
           (abs(clampAngle(START_ANGLE - angle())) < 1.0));
}

double Rot3Fold::orientationAngle(int oi, SearchStrategy ss) const {
   double delta = 0.0;
   if (ss != Mover::LOW_RES) {
      const int oh = oi;
      oi = orientation(Mover::LOW_RES);
      if (oi == 0) {   // centered on START_ANGLE
	 delta = ((oh+1)>>1) * ((oh&1) ? 1.0 : -1.0) * FINE_STEP;
      }
      else {
	 const int halfRange=int(double(ROUGH_STEP)/double(FINE_STEP) + 0.5);

	 if (oi&1) { // odd  - angles beyond START_ANGLE
	    delta = (oh * FINE_STEP) - halfRange;
	 }
	 else {      // even - angles before START_ANGLE
	    delta = halfRange - (oh * FINE_STEP);
	 }
      }
   }

   const double a = START_ANGLE + delta +
	    ((oi+1)>>1) * ((oi&1) ? 1.0 : -1.0) * ROUGH_STEP;
   return clampAngle(a);
}

std::string Rot3Fold::describeOrientation() const {
   char descrbuf[25];
   ::sprintf(descrbuf,"%s%4.0f",_grpName, _angle);
   return descrbuf;
}

#ifdef ROTPENALTY
double Rot3Fold::orientationPenalty(float pmag) const {
   const double freq = 1.5*(PI/180.0);
   // relatively small penalty
   return -pmag * 0.1 * abs(sin( freq*(angle()-START_ANGLE) ));
}
#else
double Rot3Fold::orientationPenalty(float) const {
   return 0.0;
}
#endif

std::string Rot3Fold::formatComment(std::string prefix) const
{
  char buf[200];

  ::sprintf(buf, "USER  MOD %s:%s:%s:sc=%8.3g%c  (180deg=%.3g%s)",
    prefix.c_str(),
    descr().c_str(), describeOrientation().c_str(),
    bestScore(),((bestHasBadBump())?'!':' '),
    initScore(),((initHasBadBump())?"!":""));

  return buf;
}
