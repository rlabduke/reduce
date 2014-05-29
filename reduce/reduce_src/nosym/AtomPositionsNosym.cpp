// name: AtomPositionsNosym.cpp
// author: J. Michael Word
// date written: 8/1/97
// purpose: Implementation for AtomPositions

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#if defined(_MSC_VER)
#pragma warning(disable:4786)
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
using std::cerr;
using std::endl;
using std::flush;
using std::ifstream;
using std::ios_base;

#ifdef OLD_STD_HDRS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#else
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include <cmath>
using std::sprintf;
using std::toupper;
using std::isdigit;
using std::atoi;
using std::atof;
using std::log10;
using std::pow;
using std::exp;
#endif

#include "../AtomPositions.h"
#include "DisjointSets.h"
#include "../StdResH.h"
#include "../RotMethyl.h"
#include "../RotAromMethyl.h" // for Arom methyls - Aram 08/13/12
#include "../RotDonor.h"
#include "../FlipMemo.h"
#include "../PDBrec.h"

//#include "stuff.h"
#include "../MoveableNode.h"

AtomPositions::~AtomPositions() {
  for (std::map<std::string, Mover*>::const_iterator i = _motionDesc.begin(); i != _motionDesc.end(); ++i)
    delete i->second;
  std::for_each(_excludePoints.begin(), _excludePoints.end(), DeleteObject());
  for (std::multimap<LocBlk, PDBrec*>::const_iterator it = _xyzBlocks.begin(); it != _xyzBlocks.end(); ++it)
    delete it->second;
}

void AtomPositions::put(PDBrec* r) {
  PDBrec* temp = new PDBrec(*r);
  _xyzBlocks.insert(std::make_pair(LocBlk(r->loc()), temp));
}

std::list<PDBrec*> AtomPositions::neighbors(const Point3d& p,
                                            Coord mindist, Coord maxdist) {
  return ::neighbors(p, mindist, maxdist, _xyzBlocks);
}

// ---------------------------------------------------------------
// move atom to new xyz pos.

void AtomPositions::reposition(const Point3d& prev, const PDBrec& rec) {
   LocBlk pb(prev);
   LocBlk nb(rec.loc());
   if (! (pb == nb)) { // rec has moved to a new xyz block
	   std::multimap<LocBlk, PDBrec*>::iterator it = _xyzBlocks.find(pb);
	   std::multimap<LocBlk, PDBrec*>::iterator it2;
	   PDBrec* a = NULL;
	   while (it != _xyzBlocks.end() && it->first == pb) {
		   a = it->second;
		   it2 = it;
		   ++it;
		   if (a->valid() && *a == rec) {
			   delete it2->second;
			   _xyzBlocks.erase(it2);
		   }
	   }
	   PDBrec* r = new PDBrec(rec);
	   _xyzBlocks.insert(std::make_pair(nb, r));
   }
}

// ---------------------------------------------------------------
CliqueList AtomPositions::findCliques() const {
	const float prRad = 0.0; // we just want to examine bumps here
	const int mdsz = _motionDesc.size();
	float maxVDWrad = 1.0;

	std::multimap<LocBlk, BumperPoint*> bumpbins;  //(2*mdsz);
	std::vector<MoverPtr> memo;
	memo.reserve(mdsz);
//	Vector< MoverPtr >   memo(mdsz);
	std::map<std::string, Mover*>::const_iterator it = _motionDesc.begin();
	std::string descr;

	// pack each memo into an array and add bumper points to bumpbins
	int rn=0, an=0;
	while (it != _motionDesc.end()) {
		Mover   *mx = it->second;
		//_os << "find cliques " << mx->descr() << " valid: " << mx->valid() << " complete: " << mx->isComplete() << endl;
		if (mx != NULL && mx->valid() && mx->isComplete()
			&& rn < mdsz) {
			memo.push_back(mx);
			//_os << "#res " << rn << " is " << mx->descr() << endl;
			an += mx->makebumpers(bumpbins, rn, maxVDWrad);
			++rn;
		}
		++it;
	}
	const int numMemos = rn;
//	memo.resize(numMemos);

// exclude bumper points which interact with exclude-list members

	BumperPoint* xp = NULL;
	for (std::list<BumperPoint*>::const_iterator xp_iter = _excludePoints.begin(); xp_iter != _excludePoints.end(); ++xp_iter) {
		xp = *xp_iter;
		std::list<BumperPoint*> nearxp_list = ::neighbors(xp->loc(), (Coord)0.001,
			(Coord)(xp->vdwRad() + prRad + maxVDWrad), bumpbins);
		std::list<BumperPoint*>::iterator nearxp = nearxp_list.begin();
		while (nearxp != nearxp_list.end()) {
			if ((*nearxp)->valid()
				&& ((distance2(xp->loc(), (*nearxp)->loc())
				- xp->vdwRad() - (*nearxp)->vdwRad()) < prRad)) {

				(*nearxp)->invalidate(); // turn off if bumps exclude pt
			}
			++nearxp;
		}
	}
	DisjointSets connsets(numMemos);
	std::multimap<LocBlk, BumperPoint*>::iterator bbit = bumpbins.begin();

//	add connections between memos to DisjointSet (& remember pairwise links)

	const LocBlk* key;
	const BumperPoint* bp;
	while (bbit != bumpbins.end()) {
		key = &(bbit->first);
///////////////////////////////////////////////////////////////////////////////////////
// Have trouble to use upper_bound()
//		for (; bbit != bumpbins.upper_bound(*key); ++bbit) {
		for (; bbit != bumpbins.end() && bbit->first == *key; ++bbit) {
///////////////////////////////////////////////////////////////////////////////////////
			bp = bbit->second;
			std::list<BumperPoint*> nearbp_list = ::neighbors(bp->loc(), (Coord)0.001,
				(Coord)(bp->vdwRad() + prRad + maxVDWrad), bumpbins);
			std::list<BumperPoint*>::iterator nearbp = nearbp_list.begin();
			while (nearbp != nearbp_list.end()) {
				const int ares = bp->rnum(), bres = (*nearbp)->rnum();
				if ((ares != bres) && bp->valid() && (*nearbp)->valid()
					&& (!sameres(memo[ares]->exampleAtom(), memo[bres]->exampleAtom()))
					&& ((distance2(bp->loc(), (*nearbp)->loc())
					- bp->vdwRad() - (*nearbp)->vdwRad()) < prRad)) {

#ifdef DEBUGCLIQUE
//	       _os << "-connect " << ares << " to " << bres << endl;
#endif
					connsets.connect(ares, bres); // build connected sets

					// remember who connects to whom
					if (! memo[ares]->linkExists(memo[bres]->descr())) {
						memo[ares]->addLink(memo[bres]->descr(), memo[bres]);
					}
					if (! memo[bres]->linkExists(memo[ares]->descr())) {
						memo[bres]->addLink(memo[ares]->descr(), memo[ares]);
					}
				}
				++nearbp;
			}
		}
	}

// figure out subset connectivity

   int** djss = connsets.subsets();

// package connectivity info into a CliqueList
   int i = 0, j = 0;

   CliqueList clst;
   for(j=1; j <= djss[0][0]; j++) {
      std::list<MoverPtr> c;
      for(i=1; i <= djss[j][0]; i++) {
         c.push_front(memo[ djss[j][i] ]);
      }
	  mpSeqSort(c);
      clst.insertClique(c);
   }
   freeDJsubsets(djss); // finished with array of subset indexes

// group non-interacting elements into a single list

   for(i=0; i < numMemos; i++) {
      if (connsets.singleton(i)) { clst.insertSingleton(memo[i]); }
   }
   clst.sortSingletonsByDescr();

#ifdef DEBUGCLIQUE
   _os << "Searched " << rn << " residues and examined " << an << " frontier points." << endl;
#endif

   for (std::multimap<LocBlk, BumperPoint*>::iterator iter = bumpbins.begin(); iter != bumpbins.end(); ++iter)
	   delete iter->second;

   return clst;
}


void AtomPositions::CollectBumping(const AtomDescr& ad, std::list<PDBrec*>& bumping)
{
	//std::cerr << "Collect Bumping: " << ad << " pr="  << _probeRadius << " maxVDW =" << _maxVDWFound << std::endl;
	std::list<PDBrec*> nearby = ::neighbors(ad.getAtomPos(), (Coord)0.001,
      		(Coord)(ad.getAtRadius() + _probeRadius + _maxVDWFound ), _xyzBlocks);

	PDBrec* ni = NULL;
	for(std::list<PDBrec*>::iterator it = nearby.begin(); it != nearby.end(); ++it) {
		ni = *it;
		//std:: cerr << "nearby: " << ni->getAtomDescr() << std::endl;
		double distancecutoff = (ad.getAtRadius() + ni->vdwRad() + _probeRadius);
	  if (ni->valid() && (! ni->hasProp(IGNORE))
		  && (abs(ni->occupancy()) > _occupancyCuttoff)
		  && ( distanceSquared (ad.getAtomPos(), ni->loc()) < distancecutoff * distancecutoff )) {
		  bumping.push_back(ni);
		  //cerr << "COLBUMP INSERTED: bumper = " << b.loc() << " radius = " << b.vdwRad() << " ni = " << (ni->data()).getAtomDescr() << endl;
	     }
	}
}

// ---------------------------------------------------------------
// create possible orientations for H atoms on waters (identified elsewhere)
// and store these Hs in the xyz table
void AtomPositions::generateWaterPhantomHs(std::list<PDBrec*>& waters) {
	char descrbuf[32];
	const int MaxAccDir = 25;
	struct AccDirection {
		AccDirection() : _nam(""), _gap(999.9) {}
		AccDirection(const std::string& s, const Point3d& p, float g)
			: _nam(s), _loc(p), _gap(g) {}
		std::string  _nam;
		Point3d _loc;
		float   _gap;
	} nearbyA[MaxAccDir];

	const ElementInfo& elemHOd = * ElementInfo::StdElemTbl().element("HOd");

	PDBrec* a = NULL;
	for(std::list<PDBrec*>::const_iterator it = waters.begin(); it != waters.end(); ++it) {
		a = *it;
		int i = 0, nAcc = 0;

		std::list<PDBrec*> nearby_list = neighbors(a->loc(), 0.001, 4.0);
		PDBrec* rec = NULL;
		for(std::list<PDBrec*>::const_iterator nearby = nearby_list.begin(); nearby != nearby_list.end(); ++nearby) {
			rec = *nearby;

			if (rec->hasProp(ACCEPTOR_ATOM)
				|| FlipMemo::isHBDonorOrAcceptorFlipped(*rec, _useXplorNames, _useOldNames, _bbModel)) {

				double HBoverlap = distance2(a->loc(), rec->loc())
					- ( elemHOd.explRad() + rec->vdwRad()
					+ 1.0 /*i.e. the O-H bond length*/ );

				if (HBoverlap < -0.01
					&& (abs(rec->occupancy()) > _occupancyCuttoff)
					&& interactingConfs(*a, *rec, _onlyA)) {

					// Now we make a table of all the locations of each HB
					// acceptor nearby keeping ONLY one for each aromatic ring
					// (the one with the most negative gap).

					bool foundMatchingAromAtom = FALSE;
					bool isAromRingAtom = StdResH::ResXtraInfo().atomHasAttrib(
						rec->resname(), rec->atomname(), AROMATICFLAG);

					if (!UseSEGIDasChain) {
					  ::sprintf(descrbuf, "%-3.3s%c%-3.3s%4d%c%-2.2s",
						(isAromRingAtom ? "/R/" : ""), rec->alt(),
						rec->resname(), rec->resno(),
						rec->insCode(), rec->chain());
					}
					else {
					  ::sprintf(descrbuf, "%-3.3s%c%-3.3s%4d%c%-4.4s",
						(isAromRingAtom ? "/R/" : ""), rec->alt(),
						rec->resname(), rec->resno(),
						rec->insCode(), rec->chain());
					}
					std::string resDescr = std::string(descrbuf);

					if (isAromRingAtom) {
						for (i=0; i < nAcc; i++) {
							if (resDescr == nearbyA[i]._nam) {
								if (HBoverlap < nearbyA[i]._gap) {
									nearbyA[i]._nam = resDescr;
									nearbyA[i]._loc = rec->loc();
									nearbyA[i]._gap = HBoverlap;
								}
								foundMatchingAromAtom = TRUE;
								break;
							}
						}
					}
					if (!foundMatchingAromAtom && nAcc < MaxAccDir) {
						nearbyA[nAcc]._nam = resDescr;
						nearbyA[nAcc]._loc = rec->loc();
						nearbyA[nAcc]._gap = HBoverlap;
						nAcc++;
					}
				}
			}
		}
		// each of these acceptors is a possible direction for a hydrogen

		for(i = 0; i < nAcc; i++) {
			PDBrec* pHatom = new PDBrec();
			a->clone(pHatom); // duplicate & modify Oxygen

#define BEST_HBOND_OVERLAP 0.6
			double waterOHbondLen = 1.0 +
				std::max(-1.0, std::min(0.0, nearbyA[i]._gap + BEST_HBOND_OVERLAP));

			pHatom->loc((nearbyA[i]._loc - a->loc()).scaleTo(waterOHbondLen)
				+ a->loc());
			pHatom->atomname(" H? ");

			pHatom->elem(elemHOd);

			pHatom->atomno(0);
			pHatom->occupancy(_occupancyCuttoff + 0.01);
			pHatom->tempFactor(99.99);
			pHatom->elemLabel(" H");
			pHatom->chargeLabel("  ");
			pHatom->annotation("dummy atom");

			_xyzBlocks.insert(std::make_pair(LocBlk(pHatom->loc()), pHatom));
#ifdef DUMPPHANH
			//if (waterOHbondLen < 1.0)
			_os << "{" << a->recName() << "}P "
				<< a->x()<< ", " << a->y() << ", " << a->z()
				<< endl
				<< "{" << pHatom->recName() << "}L "
				<< pHatom->x()<< ", " << pHatom->y() << ", " << pHatom->z()
				<< endl;
#endif
		}
	}
}

// ---------------------------------------------------------------
// score for a single atom a at point p based on interactions with
// nearby atoms (ignoring excluded atoms)

double AtomPositions::atomScore(const PDBrec& a, const Point3d& p,
								float nearbyRadius, const std::list<PDBrec*>& exclude,
								const DotSph& dots, float pRad,	bool onlyBumps,
								float &bumpSubScore, float &hbSubScore,	bool &hasBadBump) {

	bumpSubScore = 0.0;
	hbSubScore   = 0.0;
	hasBadBump   = FALSE;

	if ((! a.valid()) || a.hasProp(IGNORE)
		||  (! visableAltConf(a, _onlyA))
		||  (abs(a.occupancy()) <= _occupancyCuttoff) ) {
		return 0.0;
	}

	DotsForAtom * dotsToCount = 0;
	//bool output = false;

	//if (output) std::cerr << "Scoring atom: " << a.loc() << " " << a.recName() << std::endl;
	if ( scoreAtomsAndDotsInAtomsToScoreVector_ || scoreAtomsInAtomsInHighOrderOverlapList_ )
	{

		AtomDescr atomBeingScored = a.getAtomDescr();

		//for ( std::list< AtomDescr >::iterator threeWayIter = _atomsIn3WayOverlap->begin();
		//	threeWayIter != _atomsIn3WayOverlap->end(); ++threeWayIter )
		//{
		//	std::cerr << "Atom Score: atoms in 3 way overlap: " << *threeWayIter << std::endl;
		//}

		if (	scoreAtomsAndDotsInAtomsToScoreVector_ )
		{
			assert( atoms_to_score_ptr_ != 0 );
			for (int ii = 0; ii < atoms_to_score_ptr_->size(); ++ii )
			{
				//std::cerr << "Comparing against: " << (*atoms_to_score_ptr_)[ ii ].first << std::endl;
				if ( atomBeingScored == (*atoms_to_score_ptr_)[ ii ].first )
				{
					//std::cerr << "Found it!" << std::endl;
					dotsToCount = ((*atoms_to_score_ptr_)[ ii ].second );
					//for (std::vector< AtomDescr >::iterator iter = dotsToCount->begin();
					//	iter != dotsToCount->end(); ++iter)
					//{
					//	std::cerr << "count dots inside of: " << *iter <<std::endl;
					//}
					break;
				}
			}
			if (dotsToCount == 0 ){
				//std::cerr << "Did not find atom; returning 0" << std::endl;
				return 0.0;
			}
		}
		else
		{
			assert( atoms_in_high_order_overlap_ptr_ != 0 );
			if ( std::find(
				atoms_in_high_order_overlap_ptr_->begin(),
				atoms_in_high_order_overlap_ptr_->end(),
				atomBeingScored ) == atoms_in_high_order_overlap_ptr_->end() )
			{
				return 0.0;
			}
		}
	}

	//apl now that AtomPositions has decided to score this atom, collect nearby-list
	std::list<PDBrec*> nearby = this->neighbors( p, 0.001, nearbyRadius);

	std::list<PDBrec*> bumping_list; // first we collect atoms actually interacting
	PDBrec* rec = NULL;
	for (std::list<PDBrec*>::const_iterator ni = nearby.begin(); ni != nearby.end(); ++ni) {
		rec = *ni;
		if (rec->valid() && (! rec->hasProp(IGNORE))
			&& (abs(rec->occupancy()) > _occupancyCuttoff)
			&& (! (a == *rec))
			&& (interactingConfs(a, *rec, _onlyA))
			&& ( vdwGap(a, p, *rec, rec->loc()) < pRad)
			&& (! foundInList(*rec, exclude))) {
			bumping_list.push_back(rec);
		}
	}


	std::vector< PDBrec* > bumping( bumping_list.size(), NULL );
	std::copy( bumping_list.begin(), bumping_list.end(), bumping.begin() );
	//for (int ii = 0; ii < bumping.size(); ++ii )
	//{
	//	if (output) std::cerr << "bumping: " << ii << " " << bumping[ ii ] << std::endl;
	//}

	//if (false)
	//int count_atoms_to_score = bumping.size();
	//if ( scoreAtomsAndDotsInAtomsToScoreVector_ )
	//{
	//	count_atoms_to_score = 0;
	//	int first = 0;
	//	int last = bumping.size() - 1;
	//
	//	for (std::list< PDBrec*>::iterator iter = bumping_list.begin(); iter != bumping_list.end(); ++iter)
	//	{
	//		AtomDescr bumping_descr = (*iter)->getAtomDescr();
	//		if (std::find(dotsToCount->begin(), dotsToCount->end(), bumping_descr) == dotsToCount->end() )
	//		{
	//			bumping[ last ] = *iter;
	//			accept_bumping[ last ] = false;
	//			--last;
	//		}
	//		else
	//		{
	//			bumping[ first ] = *iter;
	//			++count_atoms_to_score;
	//			++first;
	//		}
	//	}
	//
	//	if ( count_atoms_to_score == 0 ) return 0.0;
	//
	//}
	//else
	//{
	//std::copy( bumping_list.begin(), bumping_list.end(), bumping.begin() );
	//}


	int HBmask = 0;
	if (a.hasProp(   DONOR_ATOM)) { HBmask |= ACCEPTOR_ATOM; }
	if (a.hasProp(ACCEPTOR_ATOM)) { HBmask |= DONOR_ATOM;    }

	const int ndots = dots.count();
	double s = 0.0;
	for (int i=0; i < ndots; i++) // then we inspect each dots interactions
	{
		if ( dotsToCount && ! dotsToCount->dotOn( i ) ) { continue; }
		const Point3d q = p + dots[i];
		const Point3d probeLoc = (pRad > 0.0)
			? p + dots[i].scaled(dots.radius()+pRad)
			: q;

		double mingap = 999.9;
		bool keepDot = FALSE;
		bool isaHB = FALSE;
		bool tooCloseHB = FALSE;
		float HBmindist = 999.9;
		PDBrec* cause = 0;

		PDBrec* b = NULL;
		//int closest_bumping = -1;
		for (int ii = 0; ii < bumping.size(); ++ii) {
			b = bumping[ ii ];
			const Point3d locb = b->loc();
			const double vdwb = b->vdwRad();

			const double squaredist = distanceSquared( probeLoc, locb );
			if ( squaredist > (vdwb + pRad ) * ( pRad + vdwb ) )
			{
				continue;
			}

			const double dist = sqrt( squaredist );
			const double probeGap = dist- pRad - vdwb;
			const double      gap = dist       - vdwb;
			if (probeGap < 0.0) {

				if (gap < mingap) {

					const bool bothCharged = a.isCharged() && b->isCharged();

					const bool chargeComplement = bothCharged
						?  ( (a.isPositive() && b->isNegative()) ||
						(a.isNegative() && b->isPositive()) )
						: FALSE;

					if( b->hasProp(HBmask) &&
						((! bothCharged) || chargeComplement)) { // H-bond
						isaHB = TRUE;

						HBmindist = bothCharged ?
					_min_charged_hb_cutoff : _min_regular_hb_cutoff;

						tooCloseHB = (gap < -HBmindist);
					}
					else { // clash or contact
						if( b->hasProp(HB_ONLY_DUMMY) ) {
							continue; // *** dummy only counts as an HBond partner
						}
						isaHB = tooCloseHB = FALSE;
					}
					cause = b;
					keepDot = TRUE;
					//closest_bumping = ii;
					mingap = gap;
				}
			}
		}

		if (keepDot) { // remove dots inside a connected atom
			PDBrec* x = NULL;
			for (std::list<PDBrec*>::const_iterator it = exclude.begin(); it != exclude.end(); ++it) {
				x = *it;
				if ( (distanceSquared(q, x->loc()) < x->vdwRad()*x->vdwRad())
					&& x->valid() && (! x->hasProp(IGNORE))
					&& interactingConfs(a, *x, _onlyA) ) {
					keepDot = FALSE;
					break;
				}
			}
		}
    	if (! keepDot ) continue;

		//if we've gotten this far, then this dot's score gets counted

    	const float GAPscale  = 0.25; // GAP score scale Factor
    	const float BUMPweight =10.0; // BUMP score weight
    	const float HBweight   = 4.0; // HBond score weight

		double dotscore = 0;
      int overlapType = 0;
      double overlap = 0.0;

      if (mingap > 0.0)      {overlap = 0.0;       overlapType =  0;}
      else if (isaHB && tooCloseHB) {
        mingap += HBmindist;
        overlap =-0.5*mingap;
        overlapType = -1;
      }
      else if (isaHB)        {overlap =-0.5*mingap;overlapType =  1;}
      else if (mingap < 0.0) {overlap =-0.5*mingap;overlapType = -1;}

      // using scale values from probe program
      if (overlapType == -1) {              // clash
        const double bmp = -BUMPweight*overlap;
        bumpSubScore += bmp;
		  dotscore = bmp;
        s += bmp;
        if (mingap <= -_bad_bump_gap_cutoff) { hasBadBump = TRUE; }
      }
      else if (overlapType == 0) {			    // contact dot
        if ((!onlyBumps)
          && ! annularDots(q, a, *cause, pRad)) {
          const double scaledGap = mingap/GAPscale;

			dotscore = scaledGap;
         dotscore = exp(-scaledGap*scaledGap);
          s += dotscore;
        }
      }
      else if (overlapType ==  1) {			    // H-bond
        if (!onlyBumps) {
          const double hbv = HBweight*overlap;
          hbSubScore += hbv;
			 dotscore = hbv;
          s += hbv;
        }
        else {   // in this case treat as a bump
          const double pseudobump = -BUMPweight*overlap;
          bumpSubScore += pseudobump;
			 dotscore = pseudobump;
          s += pseudobump;
        }
      }

  	 	//if ( scoreAtomsAndDotsInAtomsToScoreVector_ && dotscore != 0)
   	//if (output && dotscore != 0 ) cerr << "Cause " << cause->recName() << " v:" << cause->valid() << " on surface of " << a.recName() << " score: " << dotscore << " hb?" << isaHB << " tooclose?" << tooCloseHB << endl;

   }
   bumpSubScore /= dots.density();
   hbSubScore /= dots.density();
   s /= dots.density();
#ifdef DEBUGATOMSCORE
   _os << "\t:" << a.recName()
     << ": bump=" << bumpSubScore
     << ", HB=" << hbSubScore
     << ((hasBadBump==TRUE)?", BADBUMP":"")
     <<", partialScore=" << s
     << endl;
#endif
	//cerr << flush;
	//if ( _avoidAtomsIn3WayOverlap || _countAtomsIn3WayOverlapOnly ) { std::cerr << "Atom Score: " << s << std::endl; }


   return s;
}
