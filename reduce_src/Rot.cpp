// name: Rot.C
// author: J. Michael Word
// date written: 2/7/98
// purpose: Implementation for Rot

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

#include "Rot.h"
#include "AtomPositions.h"

Rot::Rot(const Point3d& a, const Point3d& b,
                     const double ang, const PDBrec& heavyAtom)
   : _p1(a), _p2(b), _heavyAtom(std::make_shared<PDBrec>(heavyAtom)), _angle(ang) {
   validateMemo();
}

double Rot::orientationAngle(int oi, SearchStrategy ss) const {
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

double Rot::determineScore(AtomPositions &xyz,
   DotSphManager& dotBucket, int nBondCutoff, float probeRadius,
   float pmag, double& penalty, float &bumpScore, float &hbScore, bool& hasBadBump) {

   bumpScore  = 0.0;
   hbScore    = 0.0;
   hasBadBump = FALSE;

   double bestScore = scoreThisAngle(xyz, dotBucket,
                         nBondCutoff, probeRadius,
                         bumpScore, hbScore,
			 hasBadBump);

   penalty = orientationPenalty(pmag);
   return bestScore;
}

double Rot::scoreThisAngle(AtomPositions &xyz, DotSphManager& dotBucket,
								 int, float probeRadius, float &bumpScore,
								 float &hbScore, bool &hasBadBump) {
	const double maxVDWrad = ElementInfo::StdElemTbl().maxExplicitRadius();

	bumpScore  = 0.0;
	hbScore    = 0.0;
	hasBadBump = FALSE;

	double scoreThisO = 0.0;
	int i = 0;
	_rot.push_front(_heavyAtom);
	for(std::list< std::shared_ptr<PDBrec> >::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
		if (!(*alst)->valid())
			continue;

		float bumpSubScore = 0.0;
		float hbSubScore   = 0.0;
		bool subBadBump    = FALSE;
		const  std::shared_ptr<PDBrec> thisAtom = *alst;

		double val = xyz.atomScore(*thisAtom, thisAtom->loc(),
      static_cast<float>(thisAtom->vdwRad() + probeRadius + maxVDWrad),
			//apl procrastinate nearby list computation until AtomPositions decides to score
			*(_bnded[i]), dotBucket.fetch(thisAtom->vdwRad()), probeRadius, FALSE,
			bumpSubScore, hbSubScore, subBadBump);

		bumpScore  += bumpSubScore;
		hbScore    += hbSubScore;
		if (subBadBump) { hasBadBump = TRUE; }

#ifdef DEBUGSUBSCORE
		cerr << "\t:" << PDBrecNAMEout(*thisAtom)
			<<":"<< describeOrientation()
			<< ": bump=" << bumpSubScore
			<< ", HB=" << hbSubScore
			<< ((subBadBump==TRUE)?", BADBUMP":"")
			<< "\t"
			//      << thisAtom.loc()
			<< endl;
#endif
		scoreThisO += val;
		i++;
	}
	_rot.pop_front();

	initializeScoreIfNotSet(scoreThisO, hasBadBump);
	return scoreThisO;
}

void Rot::setHydAngle(double newAng, AtomPositions &xyz) {
   const double oldAng = angle();
   for(std::list< std::shared_ptr<PDBrec> >::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
      setHydAngle(**hydlist, oldAng, newAng, xyz);
//cerr  << "Rot::setHydAngle[" << (*hydlist).atomname() << "] "
//      << oldAng << " deg --> "
//      << newAng << " deg" << endl;
   }
   angle(newAng);
}

void Rot::setHydAngle(PDBrec& theAtom, double oldAng, double newAng,
                        AtomPositions &xyz) {
   if (abs(newAng-oldAng) > 0.0001) {
      Point3d lastloc = theAtom.loc();
      theAtom.loc( lastloc.rotate(newAng-oldAng, _p2, _p1) );
      xyz.reposition(lastloc, theAtom);
   }
}

void Rot::dropBondedFromBumpingListForPDBrec(
	std::list< std::shared_ptr<PDBrec> > & bumping,
	std::shared_ptr<PDBrec> atom,
	int //unused nBondCutoff
) const
{
	//std::cerr << "Attempting to drop bumping for atom " << atom->getAtomDescr() << " in Rot" << std::endl;
	int atom_id = findAtom( atom );

	for (std::list< std::shared_ptr<PDBrec> >::iterator iter = bumping.begin();
		iter != bumping.end(); )
	{
		std::list< std::shared_ptr<PDBrec> >::iterator iter_next = iter;
		//std::cerr << "Comparing: " << (*iter) << " " << (*iter)->getAtomDescr() << std::endl;
		++iter_next;

		for (std::list< std::shared_ptr<PDBrec> >::const_iterator constiter = _bnded[ atom_id]->begin();
			constiter != _bnded[ atom_id]->end(); ++constiter)
		{
			if ( (*constiter)->getAtomDescr() == (*iter)->getAtomDescr() )
			{
				bumping.erase( iter );
				//std::cerr << "DROPPED!" << std::endl;
				break;
			}
		}
		iter = iter_next;
	}
}

int Rot::findAtom(std::shared_ptr<PDBrec> atom ) const
{
	if ( atom == _heavyAtom )
	{
		return 0;
	}

	int countAtomsSeen = 1;
	for (std::list< std::shared_ptr<PDBrec> >::const_iterator iter = _rot.begin(); iter != _rot.end(); ++iter)
	{
		if ( atom == *iter )
			return countAtomsSeen;
		++countAtomsSeen;
	}
	std::cerr << "Critical error in Rot::findAtom( " << atom << ").  Could not find atom. " << std::endl;
	for (std::list< std::shared_ptr<PDBrec> >::const_iterator iter = _rot.begin(); iter != _rot.end(); ++iter)
	{
		std::cerr << "_rot: " << *iter << std::endl;
	}

	exit(2);
        return 0; // to avoid warnings
}
