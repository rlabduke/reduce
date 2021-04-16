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

int Rot::makebumpers(std::multimap<LocBlk, std::shared_ptr<BumperPoint> >& bblks,
						  int rn, float& maxVDWrad) {
	int an = 0;
	std::shared_ptr<BumperPoint> bp;
	for(std::list< std::shared_ptr<PDBrec> >::const_iterator it = _rot.begin(); it != _rot.end(); ++it) {
		std::shared_ptr<PDBrec> a = *it;
		for (double theta = -scanAngle; theta < scanAngle; theta += dtheta) {
			Point3d p(a->loc().rotate(theta, _p2, _p1));
			bp = std::make_shared<BumperPoint>(p, rn, an++, a->vdwRad());
			bblks.insert(std::make_pair(LocBlk(p), bp));
			if (a->vdwRad() > maxVDWrad) { maxVDWrad = a->vdwRad(); }
		}
	}
	return an;
}

std::list<AtomDescr> Rot::getAtDescOfAllPos(float &maxVDWrad) {
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
		}
	}
	theList.sort();
	theList.unique();
	return theList;
}

bool Rot::setOrientation(int oi, float delta, AtomPositions &xyz,
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
      static_cast<float>(thisAtom->vdwRad() + maxVDWrad),
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
