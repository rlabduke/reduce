// name: RotDonor.C
// author: J. Michael Word
// date written: 2/7/98
// purpose: Implementation for RotDonor

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
#include <stdlib.h>
#else
#include <cstdio>
#include <cctype>
#include <cmath>
#include <cstdlib>
using std::exit;
#endif

#include "RotDonor.h"
#include "AtomPositions.h"

// for isHBDonorOrAcceptorFlipped
#include "FlipMemo.h"

#define START_ANGLE 180.0
#define ROUGH_STEP   10
#define FINE_STEP     1

#define OVERLAPCUTOFF_1    -0.001
#define OVERLAPCUTOFF_3    -0.7
#define OVERLAPCUTOFF_5    -0.9

RotDonor::RotDonor(const Point3d& a, const Point3d& b,
                   const double ang, const PDBrec& heavyAtom)
				   : _p1(a), _p2(b), _heavyAtom(heavyAtom), _angle(ang), _nori(0) {
	validateMemo();
}

RotDonor::~RotDonor()
{
	for (std::vector< std::list< PDBrec * > * >::iterator iter = _bnded.begin(); iter != _bnded.end(); ++iter)
	{
		std::for_each( (*iter)->begin(), (*iter)->end(), DeleteObject());
		delete (*iter);
	}
	std::for_each(_rot.begin(), _rot.end(), DeleteObject());
}


void RotDonor::finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel, 
                        AtomPositions &xyz, DotSphManager& dotBucket) {
	int i = 0;

	const double origAngle = angle();

	if (isComplete()) {

		// pre-build lists of bonded atoms

		const double approxNbondDistLimit = 3.0 + 0.5*nBondCutoff; // just a rule of thumb
		i = 0;
		_rot.push_front(&_heavyAtom);
		for (std::list<PDBrec*>::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
			const PDBrec* thisAtom = *alst;
			std::list<PDBrec*> *temp = new std::list<PDBrec*>();
			std::list<PDBrec*> neighborList = xyz.neighbors(thisAtom->loc(), thisAtom->covRad(),
				approxNbondDistLimit);
			bondedList(*thisAtom, neighborList, nBondCutoff, _rot, temp);
			_bnded.push_back(temp);
		}
		_rot.pop_front();
		// -------------------------------------------------------
		// determine all the angles we will consider

		// first identify nearby acceptors

		std::list<double> angs;

		const ElementInfo& elemHpol = * ElementInfo::StdElemTbl().element("Hpol");

		// O-H & N-H bond len == 1.0, S-H bond len == 1.3
		const double XHbondlen = (_heavyAtom.elem().atno() == 16) ? 1.3 : 1.0;

		std::list<PDBrec*> nearby_list = xyz.neighbors(_heavyAtom.loc(), 0.001, 4.0);
		PDBrec* rec = NULL;
		for(std::list<PDBrec*>::const_iterator nearby = nearby_list.begin(); nearby != nearby_list.end(); ++nearby) {
			rec = *nearby;
			if (rec->hasProp(ACCEPTOR_ATOM)
				|| FlipMemo::isHBDonorOrAcceptorFlipped(*rec, useXplorNames, useOldNames, bbModel)) {

				double HBoverlap = distance2(_heavyAtom.loc(), rec->loc())
					- ( elemHpol.explRad() + rec->vdwRad()
					+ XHbondlen );

				if (HBoverlap < OVERLAPCUTOFF_1
					&& ! diffAltLoc(_heavyAtom, *rec)) {

					const double ang = origAngle +
						dihedral(rec->loc(), _p1, _p2, (*(_rot.begin()))->loc());
//						dihedral((*nearby).loc(), _p1, _p2, (*_rot).loc());

					angs.push_front(clampAngle(ang));
					if (HBoverlap < OVERLAPCUTOFF_3) {
						angs.push_front(clampAngle(ang - 1.5*ROUGH_STEP));
						angs.push_front(clampAngle(ang + 1.5*ROUGH_STEP));
					}
					if (HBoverlap < OVERLAPCUTOFF_5) {
						angs.push_front(clampAngle(ang - 3.0*ROUGH_STEP));
						angs.push_front(clampAngle(ang + 3.0*ROUGH_STEP));
					}
				}

			}
		}

		// -------------------------------------------------------
		// next, look for least-bumping angle not near an acceptor

		double leastBumpAngle = START_ANGLE;
		double bestScore = LowestMoverScore;

		for (i=0; i < (360/ROUGH_STEP); i++) {
			const double oldTheta = angle();
			const double theta = START_ANGLE +
				((i+1)>>1) * ((i&1) ? 1.0 : -1.0) * ROUGH_STEP;
			bool tryThisAngle = TRUE;
			for(std::list<double>::iterator chi = angs.begin(); chi != angs.end(); ++chi) {
				const double delta = abs(clampAngle(clampAngle(theta) - (*chi)));
				if (delta < ROUGH_STEP) {
					tryThisAngle = FALSE;
					break;
				}
			}

			if (tryThisAngle) {
				for(std::list<PDBrec*>::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
					setHydAngle(**hydlist, oldTheta, theta, xyz);
				}
				angle(theta);
				const double bs = bumpsThisAngle(xyz, dotBucket);
				if (bs > bestScore) { bestScore = bs; leastBumpAngle = theta; }
			}
		}

		// re-position to least bump angle
		for(std::list<PDBrec*>::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
			setHydAngle(**hydlist, angle(), leastBumpAngle, xyz);
		}
		angle(leastBumpAngle);
		angs.push_front(clampAngle(leastBumpAngle));

		// -------------------------------------------------------
		// merge both sets of angles into the orientation angle array

		_nori = angs.size();
		_oang.resize(_nori);
		int j = 0;
		angs.sort();
		angs.reverse();
		for(std::list<double>::iterator chi = angs.begin(); chi != angs.end(); ++chi) {
			if (j > 0) {
				const double delta = abs(_oang[j-1] - (*chi));
				if ((delta >  (0.5*ROUGH_STEP))
					&& (*chi  >=((0.5*ROUGH_STEP)-180.0)) ) {
					_oang[j++] = *chi;
				}
			}
			else { _oang[j++] = *chi; }
		}
		_nori = j; // <-- the actual length after the merge operation
		_oang.resize(_nori);
	}
}

int RotDonor::makebumpers(std::multimap<LocBlk, BumperPoint*>& bblks,
						  int rn, float& maxVDWrad) {
	int an = 0;
	const double dtheta = 10.0; // fineness of rotation angle scan
	const double scanAngle = 180.0;
	BumperPoint* bp;
	for(std::list<PDBrec*>::const_iterator it = _rot.begin(); it != _rot.end(); ++it) {
		PDBrec* a = *it;
		for (double theta = -scanAngle; theta < scanAngle; theta += dtheta) {
			Point3d p(a->loc().rotate(theta, _p2, _p1));
			bp = new BumperPoint(p, rn, an++, a->vdwRad());
			bblks.insert(std::make_pair(LocBlk(p), bp));
			if (a->vdwRad() > maxVDWrad) { maxVDWrad = a->vdwRad(); }
		}
	}
	return an;
}

std::list<AtomDescr> RotDonor::getAtDescOfAllPos(float &maxVDWrad) {
	std::list<AtomDescr> theList;
	
	AtomDescr ad_heavy(_heavyAtom.loc(), _heavyAtom.resno(), _heavyAtom.vdwRad());
	ad_heavy.setOriginalAtomPtr( & _heavyAtom );
	theList.push_back(ad_heavy);		//ANDREW adding heavyAtom to the list of bumpers returned
	
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
		}   
	}
	theList.sort();
	theList.unique();
	return theList;
}

int RotDonor::numOrientations(SearchStrategy ss) const {
   return (ss==Mover::LOW_RES) ?
      _nori :
      (int)(2 * (double(ROUGH_STEP)/double(FINE_STEP)) + 0.5);
}

bool RotDonor::setOrientation(int oi, float delta, AtomPositions &xyz,
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

double RotDonor::orientationAngle(int oi, SearchStrategy ss) const {
   double delta = 0.0;
   if (ss != Mover::LOW_RES) {
      const int oh = oi;
      oi = orientation(Mover::LOW_RES);
      delta = ((oh+1)>>1) * ((oh&1) ? 1.0 : -1.0) * FINE_STEP;
   }

   const double a = ((oi >= 0 && oi < _nori) ? _oang[oi] : START_ANGLE)
      + delta;

   return clampAngle(a);
}

std::string RotDonor::describeOrientation() const {
   char descrbuf[25];
   ::sprintf(descrbuf, "   rot %4.0f", _angle);
   return descrbuf;
}

double RotDonor::determineScore(AtomPositions &xyz,
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

double RotDonor::scoreThisAngle(AtomPositions &xyz,	DotSphManager& dotBucket,
								int, float probeRadius,	float &bumpScore, float &hbScore,
								bool &hasBadBump) {
	const double maxVDWrad = ElementInfo::StdElemTbl().maxExplicitRadius();

	bumpScore  = 0.0;
	hbScore    = 0.0;
	hasBadBump = FALSE;

	double scoreThisO = 0.0;
	int i = 0;
	_rot.push_front(&_heavyAtom);
	for(std::list<PDBrec*>::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
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
		cerr << "\t:" << thisAtom->recName()
			<<":"<< describeOrientation()
			<< ": bump=" << bumpSubScore
			<< ", HB=" << hbSubScore
			<< ((subBadBump==TRUE)?", BADBUMP":"")
			//      << ", tot=" << val
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

double RotDonor::bumpsThisAngle(AtomPositions &xyz, DotSphManager& dotBucket) {
	const double maxVDWrad = ElementInfo::StdElemTbl().maxExplicitRadius();

	double bumpsThisA = 0.0;
	int i = 0;
	_rot.push_front(&_heavyAtom);
	for(std::list<PDBrec*>::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
		float bumpSubScore = 0.0;
		float hbSubScore   = 0.0;
		bool subBadBump    = FALSE;
		const PDBrec* thisAtom = *alst;

		double val = xyz.atomScore(*thisAtom, thisAtom->loc(), 
			thisAtom->vdwRad() + maxVDWrad, 
			//apl procrastinate nearby list computation until AtomPositions decides to score
			*(_bnded[i]), dotBucket.fetch(thisAtom->vdwRad()), 0.0, TRUE,
			bumpSubScore, hbSubScore, subBadBump);

#ifdef DEBUGSUBSCORE
		cerr <<"\tminimize clashes :"
			<< thisAtom->recName()
			<<":"<< describeOrientation()
			<<": bump="<<bumpSubScore << ((subBadBump==TRUE)?", BADBUMP":"") <<"\t"
			//      << thisAtom.loc()
			<< endl;
#endif
		bumpsThisA += val;
		i++;
	}
	_rot.pop_front();

	return bumpsThisA;
}

void RotDonor::setHydAngle(double newAng, AtomPositions &xyz) {
   const double oldAng = angle();
   for(std::list<PDBrec*>::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
      setHydAngle(**hydlist, oldAng, newAng, xyz);
//cerr  << "RotDonor::setHydAngle[" << (*hydlist).atomname() << "] "
//      << oldAng << " deg --> "
//      << newAng << " deg" << endl;
   }
   angle(newAng);
}

void RotDonor::dropBondedFromBumpingListForPDBrec(
		std::list< PDBrec * > & bumping,
		PDBrec* atom,
		int //unused nBondCutoff
) const
{
	//std::cerr << "Attempting to drop bumping for atom " << atom->getAtomDescr() << " in RotD" << std::endl;
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

int RotDonor::findAtom( PDBrec* atom ) const
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
	std::cerr << "Critical error in RotDonor::findAtom( " << atom << ").  Could not find atom. " << std::endl;
	std::cerr << "&_heavyAtom: " << &_heavyAtom << std::endl;
	for (std::list< PDBrec * >::const_iterator iter = _rot.begin(); iter != _rot.end(); ++iter)
	{
		std::cerr << "_rot: " << *iter << std::endl;
	}

	exit(2);
        return 0; // to avoid warnings
}


void RotDonor::setHydAngle(PDBrec& theAtom, double oldAng, double newAng,
                        AtomPositions &xyz) {
   if (abs(newAng-oldAng) > 0.0001) {
      Point3d lastloc = theAtom.loc();
      theAtom.loc( lastloc.rotate(newAng-oldAng, _p2, _p1) );
      xyz.reposition(lastloc, theAtom);
   }
}
