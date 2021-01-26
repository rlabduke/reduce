// name: RotMethyl.C
// author: J. Michael Word
// date written: 2/7/98
// purpose: Implementation for RotMethyl

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

#include "RotMethyl.h"
#include "AtomPositions.h"

RotMethyl::RotMethyl(const Point3d& a, const Point3d& b,
                     const double ang, const PDBrec& heavyAtom)
   : Rot(a, b, ang, heavyAtom)
{
   strcpy(_grpName, ((heavyAtom.elem().atno() == 7) ?
			"NH3+   " : "methyl "));
   validateMemo();
}

void RotMethyl::finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel,
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

int RotMethyl::makebumpers(std::multimap<LocBlk, std::shared_ptr<BumperPoint> >& bblks,
						   int rn, float& maxVDWrad) {
	int an = 0;
	const double dtheta = 10.0; // fineness of rotation angle scan
	const double scanAngle = 60.0;
	std::shared_ptr<BumperPoint> bp;
    for(std::list< std::shared_ptr<PDBrec> >::const_iterator it = _rot.begin(); it != _rot.end(); ++it) {
		std::shared_ptr<PDBrec> a = *it;
		for (double theta = -scanAngle; theta < scanAngle; theta += dtheta) {
			Point3d p(a->loc().rotate(theta, _p2, _p1));
			bp = std::make_shared<BumperPoint>(p, rn, an++, a->vdwRad());
			bblks.insert(std::make_pair(LocBlk(p), bp));
//			bblks.put(LocBlk(p), BumperPoint(p, rn, an++, (*a).vdwRad()));
			if (a->vdwRad() > maxVDWrad) { maxVDWrad = a->vdwRad(); }
		}
	}
	return an;
}

std::list<AtomDescr> RotMethyl::getAtDescOfAllPos(float &maxVDWrad) {
	std::list<AtomDescr> theList;
	AtomDescr ad_heavy(_heavyAtom->loc(), _heavyAtom->resno(), _heavyAtom->vdwRad());
	ad_heavy.setOriginalAtomPtr( _heavyAtom.get() );
	theList.push_back(ad_heavy);  // ANDREW: appending _heavyAtom to the getBumpersOfAllPos function

	std::shared_ptr<PDBrec> hyds;
	for (std::list< std::shared_ptr<PDBrec> >::const_iterator it = _rot.begin(); it != _rot.end(); ++it) {
		hyds = *it;
		Point3d initHydPoint = hyds->loc();
		for (int i = 0; i < numOrientations(Mover::LOW_RES); i++)
		{
			double theta = orientationAngle(i, Mover::LOW_RES);
			AtomDescr ad_h(initHydPoint.rotate(theta - _angle, _p2, _p1), _heavyAtom->resno(), hyds->vdwRad());
			ad_h.setOriginalAtomPtr( hyds.get() );
			theList.push_back(ad_h);
			//cerr << "TEST ROTMETHYL: " << AtomDescr(initHydPoint.rotate(theta, _p2, _p1), _heavyAtom.resno(), (hyds->data()).vdwRad()) << endl;
		}
	}
	theList.sort();
	theList.unique();
	return theList;
}

int RotMethyl::numOrientations(SearchStrategy ss) const {
   return (ss==Mover::LOW_RES) ?
      int(120.0/ROUGH_STEP + 0.5) :
      int(2.0*(double(ROUGH_STEP)/double(FINE_STEP)) + 0.5);
}

bool RotMethyl::isDefaultO(int oi, SearchStrategy ss) const {
   if (ss!=Mover::LOW_RES) { oi = orientation(Mover::LOW_RES); }
   return ((oi == 0) &&
           (abs(clampAngle(START_ANGLE - angle())) < 1.0));
}

bool RotMethyl::setOrientation(int oi, float delta, AtomPositions &xyz,
                                                     SearchStrategy ss) {

   const double oldTheta = angle();
   const double    theta = orientationAngle(oi, ss) + delta;

   if (abs(theta-oldTheta) > 0.1) {
	   for(std::list< std::shared_ptr<PDBrec> >::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
	 setHydAngle(**hydlist, oldTheta, theta, xyz);
      }
      angle(theta);
   }
   rememberOrientation(oi, ss);
   return TRUE;
}

double RotMethyl::orientationAngle(int oi, SearchStrategy ss) const {
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

std::string RotMethyl::describeOrientation() const {
   char descrbuf[25];
   ::sprintf(descrbuf,"%s%4.0f",_grpName, _angle);
   return descrbuf;
}

#ifdef ROTPENALTY
double RotMethyl::orientationPenalty(float pmag) const {
   const double freq = 1.5*(PI/180.0);
   // relatively small penalty
   return -pmag * 0.1 * abs(sin( freq*(angle()-START_ANGLE) ));
}
#else
double RotMethyl::orientationPenalty(float) const {
   return 0.0;
}
#endif

std::string RotMethyl::formatComment(std::string prefix) const
{
  char buf[200];

  ::sprintf(buf, "USER  MOD %s:%s:%s:sc=%8.3g%c  (180deg=%.3g%s)",
    prefix.c_str(),
    descr().c_str(), describeOrientation().c_str(),
    bestScore(),((bestHasBadBump())?'!':' '),
    initScore(),((initHasBadBump())?"!":""));

  return buf;
}
