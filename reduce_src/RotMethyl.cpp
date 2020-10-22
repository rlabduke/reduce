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

#define START_ANGLE 180.0
#define ROUGH_STEP   30
#define FINE_STEP     1

RotMethyl::RotMethyl(const Point3d& a, const Point3d& b,
                     const double ang, const PDBrec& heavyAtom)
   : _p1(a), _p2(b), _heavyAtom(heavyAtom), _angle(ang) {
   strcpy(_grpName, ((heavyAtom.elem().atno() == 7) ?
			"NH3+   " : "methyl "));
   validateMemo();
}

void RotMethyl::finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel,
						 AtomPositions &xyz, DotSphManager& dotBucket) {

	if (isComplete()) {

		// pre-build lists of bonded atoms

		const double approxNbondDistLimit = 3.0 + 0.5*nBondCutoff; // just a rule of thumb

		_rot.push_front(&_heavyAtom);
		for(std::list<PDBrec*>::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
			const PDBrec* thisAtom = *alst;
			if (thisAtom->valid()) {
			  std::list<PDBrec*>* temp = new std::list<PDBrec*>();
			  std::list<PDBrec*> neighborList = xyz.neighbors(thisAtom->loc(), thisAtom->covRad(),
				approxNbondDistLimit);
			  bondedList(*thisAtom, neighborList, nBondCutoff, _rot, temp);
			  _bnded.push_back(temp);
            }
		}
		_rot.pop_front();
	}
}

int RotMethyl::makebumpers(std::multimap<LocBlk, BumperPoint*>& bblks,
						   int rn, float& maxVDWrad) {
	int an = 0;
	const double dtheta = 10.0; // fineness of rotation angle scan
	const double scanAngle = 60.0;
	BumperPoint* bp;
    for(std::list<PDBrec*>::const_iterator it = _rot.begin(); it != _rot.end(); ++it) {
		PDBrec* a = *it;
		for (double theta = -scanAngle; theta < scanAngle; theta += dtheta) {
			Point3d p(a->loc().rotate(theta, _p2, _p1));
			bp = new BumperPoint(p, rn, an++, a->vdwRad());
			bblks.insert(std::make_pair(LocBlk(p), bp));
//			bblks.put(LocBlk(p), BumperPoint(p, rn, an++, (*a).vdwRad()));
			if (a->vdwRad() > maxVDWrad) { maxVDWrad = a->vdwRad(); }
		}
	}
	return an;
}

std::list<AtomDescr> RotMethyl::getAtDescOfAllPos(float &maxVDWrad) {
	std::list<AtomDescr> theList;
	AtomDescr ad_heavy(_heavyAtom.loc(), _heavyAtom.resno(), _heavyAtom.vdwRad());
	ad_heavy.setOriginalAtomPtr( & _heavyAtom );
	theList.push_back(ad_heavy);  // ANDREW: appending _heavyAtom to the getBumpersOfAllPos function

	PDBrec* hyds = NULL;
	for (std::list<PDBrec*>::const_iterator it = _rot.begin(); it != _rot.end(); ++it) {
		hyds = *it;
		Point3d initHydPoint = hyds->loc();
		for (int i = 0; i < numOrientations(Mover::LOW_RES); i++)
		{
			double theta = orientationAngle(i, Mover::LOW_RES);
			AtomDescr ad_h(initHydPoint.rotate(theta - _angle, _p2, _p1), _heavyAtom.resno(), hyds->vdwRad());
			ad_h.setOriginalAtomPtr( hyds );
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
	   for(std::list<PDBrec*>::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
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

double RotMethyl::determineScore(AtomPositions &xyz,
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

double RotMethyl::scoreThisAngle(AtomPositions &xyz, DotSphManager& dotBucket,
								 int, float probeRadius, float &bumpScore,
								 float &hbScore, bool &hasBadBump) {
	const double maxVDWrad = ElementInfo::StdElemTbl().maxExplicitRadius();

	bumpScore  = 0.0;
	hbScore    = 0.0;
	hasBadBump = FALSE;

	double scoreThisO = 0.0;
	int i = 0;
	_rot.push_front(&_heavyAtom);
	for(std::list<PDBrec*>::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
		if (!(*alst)->valid())
			continue;

		float bumpSubScore = 0.0;
		float hbSubScore   = 0.0;
		bool subBadBump    = FALSE;
		const PDBrec* thisAtom = *alst;

		double val = xyz.atomScore(*thisAtom, thisAtom->loc(),
			thisAtom->vdwRad() + probeRadius + maxVDWrad,
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

void RotMethyl::setHydAngle(double newAng, AtomPositions &xyz) {
   const double oldAng = angle();
   for(std::list<PDBrec*>::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
      setHydAngle(**hydlist, oldAng, newAng, xyz);
//cerr  << "RotMethyl::setHydAngle[" << (*hydlist).atomname() << "] "
//      << oldAng << " deg --> "
//      << newAng << " deg" << endl;
   }
   angle(newAng);
}

void RotMethyl::setHydAngle(PDBrec& theAtom, double oldAng, double newAng,
                        AtomPositions &xyz) {
   if (abs(newAng-oldAng) > 0.0001) {
      Point3d lastloc = theAtom.loc();
      theAtom.loc( lastloc.rotate(newAng-oldAng, _p2, _p1) );
      xyz.reposition(lastloc, theAtom);
   }
}

void RotMethyl::dropBondedFromBumpingListForPDBrec(
	std::list< PDBrec * > & bumping,
	PDBrec* atom,
	int //unused nBondCutoff
) const
{
	/*
	int atom_id = findAtom( atom );
	for (std::list< PDBrec* >::iterator iter = bumping.begin();
		iter != bumping.end(); )
	{
		std::list< PDBrec* >::iterator iter_next = iter;
		++iter_next;
		if ( std::find( _bnded[ atom_id]->begin(), _bnded[ atom_id]->end(), *iter ) != _bnded[ atom_id]->end() )
		{
			bumping.erase( iter );
		}
		iter = iter_next;
	}
	*/
	//std::cerr << "Attempting to drop bumping for atom " << atom->getAtomDescr() << " in RotMethyl" << std::endl;
	int atom_id = findAtom( atom );
	//for (std::list< PDBrec* >::const_iterator iter = _bnded[ atom_id]->begin();
	//	iter != _bnded[ atom_id]->end(); ++iter)
	//{
	//	std::cerr << "Bonded to atom: " << (*iter) << " " << (*iter)->getAtomDescr()  << std::endl;
	//}

	for (std::list< PDBrec* >::iterator iter = bumping.begin();
		iter != bumping.end(); )
	{
		std::list< PDBrec* >::iterator iter_next = iter;
		//std::cerr << "Comparing: " << (*iter) << " " << (*iter)->getAtomDescr() << std::endl;
		++iter_next;

		for (std::list< PDBrec* >::const_iterator constiter = _bnded[ atom_id]->begin();
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


int RotMethyl::findAtom( PDBrec* atom ) const
{
	if ( atom == & _heavyAtom )
	{
		return 0;
	}

	int countAtomsSeen = 1;
	for (std::list< PDBrec * >::const_iterator iter = _rot.begin(); iter != _rot.end(); ++iter)
	{
		if ( atom == *iter )
			return countAtomsSeen;
		++countAtomsSeen;
	}
	std::cerr << "Critical error in RotMethyl::findAtom( " << atom << ").  Could not find atom. " << std::endl;
	std::cerr << "&_heavyAtom: " << &_heavyAtom << std::endl;
	for (std::list< PDBrec * >::const_iterator iter = _rot.begin(); iter != _rot.end(); ++iter)
	{
		std::cerr << "_rot: " << *iter << std::endl;
	}

	exit(2);
        return 0; // to avoid warnings
}
