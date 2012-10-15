// name: AtomPositions.C
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

#include "AtomPositions.h"
#include "DisjointSets.h"
#include "StdResH.h"
#include "RotMethyl.h"
#include "RotDonor.h"
#include "FlipMemo.h"
#include "PDBrec.h"

//#include "stuff.h"
#include "MoveableNode.h"

int AtomPositions::forceOrientations(const std::string& ofilename, std::vector<std::string>& orientNotes) {
	int numInNonDefaultOrientations = 0;
	ifstream is(ofilename.c_str());
	if (! is) {
		cerr << "ERROR in Fix: can't open orientation file \""
			<< ofilename << "\"." << endl;
		return 0;
	}

	const int buflen = 200;
	char buf[buflen+1];

	while (is.getline(buf, buflen)) {
		char orient = ' ';
		float angle = 0.0;
		int oval    = 0;
		char *p = buf;
		for (; *p && (p[0] != ':') && (p[0] != '#'); p++) {
			if (orient == ' ') {
				if (toupper(*p) == 'O') { orient = 'O'; }
				if (toupper(*p) == 'F') { orient = 'F'; }
				if (toupper(*p) == 'R') { orient = 'R'; }
			}
			else if (orient == 'F') {
				if (*p == '+' || *p == '-' || isdigit(*p)) {
					oval = atoi(p);
					while (*p && (p[0] != ':')) { p++; }
					break;
				}
			}
			else if (orient == 'R') {
				if (*p == '+' || *p == '-' || isdigit(*p)) {
					angle = atof(p);
					while (*p && (p[0] != ':')) { p++; }
					break;
				}
			}
		}
		std::string descr;
		if (p[0] == ':') {
			char *q = ++p;
			while (*p && (p[0] != ':')) { p++; }
			p[0] = '\0';
			descr = q;
		}
		Mover *mx;
		if (descr.empty()) { continue; } // incomplete line or comment
		else {
			std::map<std::string, Mover*>::const_iterator iter = _motionDesc.find(descr);
			if (iter != _motionDesc.end())
				mx = iter->second;
			else
				mx = NULL;
		}
		if (mx != NULL && mx->valid() /* && mx->isComplete() */ ) {
			if (orient == 'O' && (mx->canRotate()
				|| (mx->canFlip() && (mx->numOrientations() <= 2)) )) {
				if (_outputNotice) {
					cerr << "FIX orientation=ORIG:" <<descr<< ":" << endl;
				}
				::sprintf(buf, ( "USER  MOD Fix    :%s:%s: (o)" ),
					mx->descr().c_str(), mx->describeOrientation().c_str());
				orientNotes.push_back(buf);
				mx->makeNonAdjustable();
			}
			else if (orient == 'O' && mx->canFlip()) {
				mx->limitOrientations(FALSE);
				if (_outputNotice) {
					cerr << "LIM orientation=ORIG:" <<descr<< ": "
						<< mx->numOrientations() << " states" << endl;
				}
				::sprintf(buf, ( "USER  MOD Limit  :%s:           : (o, %d states)" ),
					mx->descr().c_str(),
					mx->numOrientations());
				orientNotes.push_back(buf);
			}
			else if (orient == 'F'
				&& mx->canFlip() && (mx->numOrientations() <= 2)) {

				mx->setOrientation(1, *this);

				if (_outputNotice) {
					cerr << "FIX orientation=FLIP:" <<descr<< ":" << endl;
				}
				::sprintf(buf, ( "USER  MOD Fix    :%s:%s: (f)" ),
					mx->descr().c_str(), mx->describeOrientation().c_str());
				orientNotes.push_back(buf);
				mx->makeNonAdjustable();
				numInNonDefaultOrientations++;
			}
			else if (orient == 'F' && mx->canFlip()) {
				if ((oval > 0) && (oval <= mx->numOrientations())) {
					mx->setOrientation(oval-1, *this);

					if (_outputNotice) {
						cerr << "FIX orientation=  "<<oval<<" :" <<descr<< ":" << endl;
					}
					::sprintf(buf, ( "USER  MOD Fix    :%s:%s: (f %d)" ),
						mx->descr().c_str(),
						mx->describeOrientation().c_str(), oval);
					orientNotes.push_back(buf);
					mx->makeNonAdjustable();
					numInNonDefaultOrientations++;
				}
				else if (oval == 0) {
					mx->limitOrientations(TRUE);

					if (_outputNotice) {
						cerr << "LIM orientation=FLIP:" <<descr<< ": "
							<< mx->numOrientations() << " states" << endl;
					}
					::sprintf(buf, ( "USER  MOD Limit  :%s:FLIP       : (f, %d states)" ),
						mx->descr().c_str(),
						mx->numOrientations());
					orientNotes.push_back(buf);
				}
				else {
					cerr << "ERROR in Fix:" << descr
						<< ": can't set orientation to "
						<< oval << ", range = 1.."
						<< mx->numOrientations() << endl;
				}
			}
			else if (orient == 'R' && mx->canRotate()) {
				mx->setHydAngle(clampAngle(angle), *this);
				if (_outputNotice) {
					cerr << "FIX rot ang="<<clampAngle(angle)
						<< " deg :" <<descr<< ":" << endl;
				}
				::sprintf(buf, ( "USER  MOD Fix    :%s:%s: (r %.0f)" ),
					mx->descr().c_str(),
					mx->describeOrientation().c_str(), angle);
				orientNotes.push_back(buf);
				mx->makeNonAdjustable();
				numInNonDefaultOrientations++;
			}
			else {
				cerr << "ERROR in Fix:" << descr << ": don't know how to '"
					<< buf << "'" << endl;
			}
		}
		else { cerr << "ERROR in Fix:" << descr
			<< ": not a valid adjustable group." << endl; 
		}
	}
	return numInNonDefaultOrientations;
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
void AtomPositions::insertRot(const PDBrec& hr,
                const PDBrec& c1, const PDBrec& c2, const PDBrec& c3,
		bool doOHSH, bool doNH3, bool doMethyl) {
	
	if ((! visableAltConf(hr, _onlyA))
		|| (! visableAltConf(c1, _onlyA))
		|| (! visableAltConf(c2, _onlyA))) { return; }

	char descrbuf[30];

	::sprintf(descrbuf, "%-2.2s%4s%c%-3.3s%-4.4s%c",
		hr.chain(), hr.Hy36resno(), hr.insCode(),
		hr.resname(), c1.atomname(), hr.alt());
	std::string descr = descrbuf;

	const double dang = dihedral(hr.loc(), c1.loc(), c2.loc(), c3.loc());

// cerr << "AtomPositions::insertRot(dang=" << dang << ") "
//      << hr.chain() << " " <<  hr.resno() << hr.insCode()
//      << " " << hr.resname() << " " << hr.atomname() <<  hr.alt() << endl;

	std::map<std::string, Mover*>::const_iterator iter = _motionDesc.find(descr);
	Mover *m;
	if (iter != _motionDesc.end())
		m = iter->second;
	else
		m = NULL;

	if (m == NULL) {
		if (c1.elem().atno() == 6) {
			if (doMethyl) {
				m = new RotMethyl(c1.loc(), c2.loc(), dang, c1);
				_motionDesc.insert(std::make_pair(descr, m));
			}
		}
		else if (c1.elem().atno() == 7) {
			if (doNH3) {
				m = new RotMethyl(c1.loc(), c2.loc(), dang, c1);
				_motionDesc.insert(std::make_pair(descr, m));
			}
		}
		else if (c1.hasProp(ACCEPTOR_ATOM)) {
			if (doOHSH) {
				m = new RotDonor(c1.loc(), c2.loc(), dang, c1);
				_motionDesc.insert(std::make_pair(descr, m));
			}
		}
		else {
			cerr<<"*error* insertRot(heavy atom("<<
				c1.recName()<<") not Carbon and not ACCEPTOR)"<<endl;
		}
	}
	if (m != NULL) {
		if (m->type() == Mover::ROTATE_METHYL){
			if (doMethyl || doNH3) {
				((RotMethyl*)m)->insertHatom(hr);
			}
		}
		else if (m->type() == Mover::ROTATE_DONOR){
			if (doOHSH) {
				((RotDonor*)m)->insertHatom(hr);
			}
		}
		else {
			cerr<<"*error* insertRot(hr, "<< m->type() <<")"<<endl;
		}
	}
}

// ---------------------------------------------------------------
// Make the group non-adjustable (fixed)

void AtomPositions::doNotAdjust(const PDBrec& a) {
	char descrbuf[30];

	::sprintf(descrbuf, "%-2.2s%4d%c%-3.3s%-4.4s%c",
		a.chain(), a.resno(), a.insCode(),
		a.resname(), "", a.alt());
	const std::string descr = descrbuf;

	std::map<std::string, Mover*>::const_iterator iter = _motionDesc.find(descr);
	Mover *m = NULL;
	if (iter != _motionDesc.end())
		m = iter->second;

	if (m != NULL) {
		m->makeNonAdjustable();
	}
}

// ---------------------------------------------------------------
// Insert the atoms in a residue block into one or more flip memos
// and return the alt codes for the memos (empty if res not flipped)

std::list<char> AtomPositions::insertFlip(const ResBlk& rblk) {
	std::list<char> resAlts;
	FlipMemo::altCodes(rblk, _useXplorNames, _useOldNames, _bbModel, resAlts);
	std::multimap<std::string, PDBrec*> pdb = rblk.atomIt();
	std::string key, descriptor;
	char descrbuf[30];

	for(std::list<char>::iterator alts = resAlts.begin(); alts != resAlts.end(); ++alts) {
		std::multimap<std::string, PDBrec*>::const_iterator pdbit = pdb.begin();
		PDBrec* atsq = NULL;
		while(pdbit != pdb.end()) {
			key = pdbit->first;
			if (key[0] == '-' || key[0] == '+') {++pdbit; continue; } // skip connectors

			for (; pdbit != pdb.end() && pdbit->first == key; ++pdbit) {
				atsq = pdbit->second;
				if ( visableAltConf(*atsq, _onlyA)
					&&  (atsq->alt() == ' ' || atsq->alt() == *alts) ) {
  
					::sprintf(descrbuf, "%-2.2s%4d%c%-3.3s%-4.4s%c",
						atsq->chain(), atsq->resno(), atsq->insCode(),
						atsq->resname(), "", (*alts));
					descriptor = descrbuf;

					std::map<std::string, Mover*>::const_iterator iter = _motionDesc.find(descriptor);
					Mover *m = NULL;
					if (iter != _motionDesc.end())
						m = iter->second;

					if (m == NULL) {
					   
						m = new FlipMemo(atsq->resname(), _useXplorNames, _useOldNames, _bbModel);
						m->descr( descriptor );
						_motionDesc.insert(std::make_pair(descriptor, m));
						//std::cerr << "FlipMemo constructed: " << descriptor << " " << m << std::endl;
					}
					else if (m->type() != Mover::FLIP){
						cerr<<"*error* insertFlip(rblk, "<< m->type() <<")"<<endl;
						return std::list<char>();
					}
					((FlipMemo*)m)->insertAtom(atsq);
					 
				}
			}
		}
	}
	return resAlts;
}

// Insert an atom into one or more flip memos based on the residues alt codes

void AtomPositions::insertFlip(PDBrec* hr, std::list<char> alts_list) {
	char descrbuf[30];

	if (! visableAltConf(*hr, _onlyA)) { return; }

	std::list<char>::iterator alts = alts_list.begin();
	while(alts != alts_list.end()) {
		if (hr->alt() == ' ' || hr->alt() == *alts) {
			::sprintf(descrbuf, "%-2.2s%4d%c%-3.3s%-4.4s%c",
				hr->chain(), hr->resno(), hr->insCode(),
				hr->resname(), "", (*alts));
			std::string descr = descrbuf;

			std::map<std::string, Mover*>::const_iterator iter = _motionDesc.find(descr);
			Mover *m = NULL;
			if (iter != _motionDesc.end())
				m = iter->second;

			if (m == NULL) {
				m = new FlipMemo(hr->resname(), _useXplorNames, _useOldNames, _bbModel);
				_motionDesc.insert(std::make_pair(descr, m));
			}
			else if (m->type() != Mover::FLIP){
				cerr<<"*error* insertFlip(hr, "<< m->type() <<")"<<endl;
				return;
			}
			((FlipMemo*)m)->insertAtom(hr);
		}
		++alts;
	}
}
// ---------------------------------------------------------------
void AtomPositions::finalizeMovers() {
	std::map<std::string, Mover*>::const_iterator it = _motionDesc.begin();
	int rn=0;
	while (it != _motionDesc.end()) {
		Mover *m = it->second;
		if (m != NULL) {
			
			m->descr(it->first);
			//_os << "FinalizeMovers: #res " << rn << " is " << m->descr() << endl;
			m->finalize(_nBondCutoff, _useXplorNames, _useOldNames, _bbModel, *this, _dotBucket);
			//_os << "FinalizeMovers: #res " << rn << " is " << m->descr() << endl;
		}
		++rn;
		++it;
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

// ---------------------------------------------------------------
int AtomPositions::orientSingles(const std::list<MoverPtr>& singles) {
#ifndef OLD_CXX_DEFNS
//	const long lff = _os.setf(     ios::fixed,      ios::floatfield); //2.22b ill-advised modification reverted
	ios_base::fmtflags lff = _os.setf(ios_base::fixed, ios_base::floatfield);
#else
	const long lff = _os.setf(     ios::fixed,      ios::floatfield);
#endif
#ifndef OLD_CXX_DEFNS
//	const long laf = _os.setf(     ios::right,      ios::adjustfield); //2.22b ill-advised modification reverted
	ios_base::fmtflags laf = _os.setf(ios_base::right, ios_base::adjustfield);
#else
	const long laf = _os.setf(     ios::right,      ios::adjustfield);
#endif
	const int iw = _os.width(8);
	const int ip = _os.precision(3);

	int i;

	double penalty = 0.0;
	int numInNonDefaultOrientations = 0;

	float bumpTrialScore = 0.0;
	float hbTrialScore   = 0.0;
	bool trialBadBump    = FALSE;

//	while (singles) {
	for (std::list<MoverPtr>::const_iterator iter = singles.begin(); iter != singles.end(); ++iter) {
//		MoverPtr m = *singles++;
		MoverPtr m = *iter;
		if (m && m->valid()) {
			double prevPenalty = 0.0;
			float theBumpScore = 0.0;
			float theHBScore   = 0.0;
			bool hasBadBump    = FALSE;
			const int numLowResO = m->numOrientations();
			m->setBestOrientation(0, LowestMoverScore, hasBadBump);

			// try each orientation
			for (i = 0; i < numLowResO; i++) {
				if (m->setOrientation(i, *this)) {
					const double val = m->determineScore(*this,
						_dotBucket, _nBondCutoff, _probeRadius,
						_pmag, penalty, bumpTrialScore,
						hbTrialScore, trialBadBump);
					if (val+penalty > m->bestScore()+prevPenalty) {
						m->setBestOrientation(i, val, trialBadBump);
						theBumpScore = bumpTrialScore;
						theHBScore   = hbTrialScore;
						hasBadBump   = trialBadBump;
						prevPenalty  = penalty;
					}
					if (_outputNotice && _showOrientScore) {
						_os << m->descr()
							<< "[" << i+1 << ":" << m->describeOrientation()
							<<"] bump=" << bumpTrialScore
							<< ", HB=" << hbTrialScore
							<< ", penalty=" << penalty
							<< ", bestScore=" << m->bestScore()
							<< ((trialBadBump==TRUE)?", BADBUMP":"")
							<< endl;
					}
					if (m->canFlip()) {
						m->trackFlipMaxScore(m->flipState(), val, trialBadBump);
					}
				}
				else { break; } // could not set orientation (we are done)
			}
			
			// remember and position to the best
			const int best = m->bestOrientation();
			if (m->orientation() != best) {
				m->setOrientation(best, *this);
			}
			if (m->hasHires()) {
				const int numHiResO = m->numOrientations(Mover::HIGH_RES);
				// try each high-res orientation
				for (i = 0; i < numHiResO; i++) {
					if (m->setOrientation(i, *this, Mover::HIGH_RES)) {
						const double val = m->determineScore(*this,
							_dotBucket, _nBondCutoff, _probeRadius,
							_pmag, penalty, bumpTrialScore,
							hbTrialScore, trialBadBump);

						if ((i == 0) || (val+penalty > m->bestScore()+prevPenalty)) {
							m->setBestOrientation(i, val, trialBadBump, Mover::HIGH_RES);
							theBumpScore = bumpTrialScore;
							theHBScore   = hbTrialScore;
							hasBadBump   = trialBadBump;
							prevPenalty  = penalty;
						}
						if (_outputNotice && _showOrientScore) {
							_os << ">" << m->descr()
								<< "[" << best+1 << "." << i+1 << ":"
								<< m->describeOrientation()
								<<"] bump=" << bumpTrialScore
								<< ", HB=" << hbTrialScore
								<< ", penalty=" << penalty
								<< ", bestScore=" << m->bestScore()
								<< ((trialBadBump==TRUE)?", BADBUMP":"")
								<< endl;
						}
						if (m->canFlip()) {
							m->trackFlipMaxScore(m->flipState(), val, trialBadBump);
						}
					}
					else { break; } // could not set orientation (we are done)
				}
				
				// remember and position to the best
				const int bestHires = m->bestOrientation(Mover::HIGH_RES);
				if (m->orientation(Mover::HIGH_RES) != bestHires) {
					m->setOrientation(bestHires, *this, Mover::HIGH_RES);
				}
				if (! m->isDefaultO(bestHires, Mover::HIGH_RES)) {
					numInNonDefaultOrientations++;
				}
			}
			else {
				if (! m->isDefaultO(best)) {
					numInNonDefaultOrientations++;
				}
			}
			if (m->canFlip() && (m->flipState() != 0)) {
				m->markFlipAtoms();
			}
			if (_outputNotice) {
				_os << " orientation " << best+1
					<< ":" << m->descr()
					<< ":" << m->describeOrientation()
					<< ": bump=" << theBumpScore
					<< ", HB=" << theHBScore
					<< ", total=" << m->bestScore()
					<< ((hasBadBump==TRUE)?", BADBUMP":"")
					<< endl;
			}
      }
   }
#ifndef OLD_CXX_DEFNS
//   _os.setf(lff, ios::floatfield); //2.22b ill-advised modification reverted
//   _os.setf(laf, ios::adjustfield); //2.22b ill-advised modification reverted
   _os.setf(lff, ios_base::floatfield);
   _os.setf(laf, ios_base::adjustfield);
#else
   _os.setf(lff, ios::floatfield);
   _os.setf(laf, ios::adjustfield);
#endif
   _os.width(iw);
   _os.precision(ip);
   return numInNonDefaultOrientations;
}

// ---------------------------------------------------------------
int AtomPositions::orientClique(const std::list<MoverPtr>& clique, int limit) {
#ifndef OLD_CXX_DEFNS
//	const long lff = _os.setf(     ios::fixed,      ios::floatfield); //2.22b ill-advised modification reverted
	ios_base::fmtflags lff = _os.setf(ios_base::fixed, ios_base::floatfield);
#else
	const long lff = _os.setf(     ios::fixed,      ios::floatfield);
#endif
#ifndef OLD_CXX_DEFNS
//	const long laf = _os.setf(     ios::right,      ios::adjustfield); //2.22b ill-advised modification reverted
	ios_base::fmtflags laf = _os.setf(ios_base::right, ios_base::adjustfield);
#else
	const long laf = _os.setf(     ios::right,      ios::adjustfield);
#endif
	const int iw = _os.width(8);
	const int ip = _os.precision(3); //ANDREW: 3 originally - sometimes 6 is useful to output more precision
	
	double log10ncombo = 0.0;
	int i = 0, rc = 0;
	
	if (_outputNotice) { _os << "\n Processing set"; }
	for(std::list<MoverPtr>::const_iterator s = clique.begin(); s != clique.end(); ++s) {
		const int no = (*s)->numOrientations();
		log10ncombo += log10(double(no));
		if (_outputNotice && ++i > 3) { i = 0; _os << endl<< "  "; }
		if (_outputNotice) { _os << ":" << (*s)->descr() << "[" << no << "]"; }
	}	
	const double p10 = floor(log10ncombo);
	const float 	mantis = pow(double(10.0), log10ncombo - p10);
	const int 		poweroften = int(p10);
	if(_outputNotice){
		if (log10ncombo < log10(double(1000000))) { //smallish numbers as integers
			const long ncombo = (long)floor(pow(double(10.0), log10ncombo)+0.5);
			_os << std::endl << " permutations: " << ncombo << endl;
		}
		else {
			_os << std::endl << " permutations: " << mantis <<"E" << poweroften << endl;
		}
	}
	
	rc = SearchClique(clique, limit);   	

#ifndef OLD_CXX_DEFNS
//	_os.setf(lff, ios::floatfield); //2.22b ill-advised modification reverted
//	_os.setf(laf, ios::adjustfield); //2.22b ill-advised modification reverted
	_os.setf(lff, ios_base::floatfield);
	_os.setf(laf, ios_base::adjustfield);
#else
	_os.setf(lff, ios::floatfield);
	_os.setf(laf, ios::adjustfield);
#endif
	_os.width(iw);
	_os.precision(ip);
	
	return rc;
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

int AtomPositions::SearchClique(std::list<MoverPtr> clique, int limit)
{
	const int numItems = clique.size();
	std::vector<MoverPtr>    item (numItems);
	std::vector<double> currScore;
	currScore.reserve(numItems);
	std::vector<float>   currBump;
	currBump.reserve(numItems);
	std::vector<float>  currHbond;
	currHbond.reserve(numItems);
	std::vector<bool> currBadBump;
	currBadBump.reserve(numItems);
	std::vector<std::list<AtomDescr> >   atomsInCliq( numItems );
	std::vector<std::list<AtomDescr>* >   bmpingAtoms;
	std::vector<std::list<AtomDescr > > bumpingNonBonded( numItems );
	std::vector< std::vector< float > > penalties( numItems );
	std::list< float > allPenalties;

	std::vector< std::list< AtomDescr > > atomsIn3WayOverlapForNode( numItems );
	std::vector< bool > nodeHas3WayOverlap( numItems, false);
	
	//We calculate edge scores first - unambiguous atoms that
	//are capable of simultaneously interacting with two Moveables
	//need to be kept track of for each of the Moveables so they can
	//be later excluded from the vertex scoring.
	double scoreThisO       = 0.0;
	double penalty          = 0.0;
	double worstPenalty     = 0.0;
	//  double prevWorstPenalty = 0.0;
	int i = 0, best = 0, j=0;
   _maxVDWFound = ElementInfo::StdElemTbl().maxExplicitRadius();

	// ---------------------------------------------------------------
	// initialize

	if ( ! initializeCliqueMovers( clique, item, numItems ) ) return 0;
	std::vector< int > num_states( numItems, 0 ); 
	
	//if ( _outputNotice )
	//{
	//	for (i=0; i<numItems;i++){
	//		std::cerr << "ITEM DESCR: " << item[i]->descr();
	//	}
	//}

   //int  countCliqAtoms = 0;
   
	setNumStatesForNodes( item, numItems, num_states, penalties );

	GraphToHoldScores gths( this, & _dotBucket, num_states, item ); 
	//all low-order (4-way overlap and lower) scoring is done in constructor

	gths.getPenalties( penalties, allPenalties );
	
	int numFlipableNodes = 0;
	for (int ii = 0; ii < numItems; ++ii)
	{
		if ( item[ii]->canFlip()) {++numFlipableNodes;}
	}

	//find out how many possible penalty values there are
	allPenalties.sort();
	allPenalties.unique();
	int numPenaltyValues = allPenalties.size();

 	//---------- Optimize Based on the Interaction Graph ------------// 
 	std::vector< NodeState > optimal_state_enabled( numItems ); //enabled enumeration of states
 	std::vector< int > optimal_state_original( numItems );//original enumeration of states
 	std::vector< NodeState > optimal_state( numItems ); //after determining penalty
 	float best_score_wo_penalty = -1;
 	float best_score_including_penalty = -1;
 	
	
	bool abandonedOptimization = false;
	bool firstOptimization = true;

	int numOptimizationProblemsToSolve = numPenaltyValues + 	numFlipableNodes;
	
	if (_outputNotice) 
	{
		std::cerr << " Num optimizations problems to be solved for this clique: ";
		std::cerr << numOptimizationProblemsToSolve << std::endl;
	}
	double timeLimit = limit / numOptimizationProblemsToSolve;
	
	allPenalties.reverse();
	for ( std::list< float >::const_iterator penalty_iterator = allPenalties.begin();
		penalty_iterator != allPenalties.end(); ++penalty_iterator)
	{
		//disable all states with penalties smaller (more negative) than *penalty_iterator
		gths.setAllStatesEnabled();
		for (int ii = 0; ii < numItems; ++ii)
		{
			int const ii_num_states = gths.getNumStatesForNode( ii );
			for (int jj = 0; jj < ii_num_states; ++jj)
			{
				if ( penalties[ ii ][ jj ] < *penalty_iterator )
				{
					//std::cerr << "Disabling: " << ii << " " << jj << " " << penalties[ ii ][ jj ] << " " << *penalty_iterator << std::endl;
					gths.disableStateOnNode( ii, jj );
					
				}
			}
		}
		gths.setStateDisablingCompleteForNow();
		//std::cerr << "state disabling complete" << std::endl;
		
		if ( gths.anyNodeWithAllStatesDisabled() )
		{
			std::cerr << " Skipping iteration; node with all states disabled" << std::endl;
			continue;
		}
		
		NodeAndEdgeManager* theNaEManager = NodeAndEdgeManager::getInstance();
		theNaEManager->InitializeNetwork(gths);
		if ( firstOptimization )
		{
			theNaEManager->setTimeLimit( timeLimit );
		}
		abandonedOptimization = theNaEManager->computeOptimalNetworkConfiguration();
		//theNaEManager->bruteForceOriginalGraph();
		
		if ( abandonedOptimization ) {
			if (_outputNotice )
			{
				std::cerr << " Abandoned Optimization.  Projected that " << numOptimizationProblemsToSolve;
				std::cerr << " optimizations would take more than " << limit << " seconds." << std::endl;
				std::cerr << "Increase time LIMIT to optimize this clique" << std::endl;
			}
			theNaEManager->clear();
			return -1;
		}
		
		optimal_state_enabled = theNaEManager->getOptimalNetworkState();
		for (int ii = 0; ii < numItems; ++ii)
		{
			optimal_state_original[ ii ] = 
				gths.convertEnabledState2OriginalStateEnumerationOnNode(
				ii, optimal_state_enabled[ ii ] );
			//std::cerr << "optimal enabled state for " << ii << " : " << optimal_state_enabled[ ii ] << " --> original : " << optimal_state_original[ ii ] << std::endl;
		}
		
		//std::cerr << "Score According to gths: " << gths.getScoreForStateAssignment( optimal_state_original ) << std::endl;
		
		float actualPenalty = 0;
		for (int ii = 0; ii < numItems; ++ii)
		{
			if ( ii == 0 || actualPenalty > penalties[ ii ][ optimal_state_original[ ii ]] )
			{
				actualPenalty = penalties[ ii ][ optimal_state_original[ ii ]];
			}
		}
		float score_without_penalty = theNaEManager->getScoreOfOptimalNetworkConfiguration();
		float score_including_penalty = score_without_penalty + actualPenalty;
		
		//std::cerr << "DP Iteration: penalty threshold " << *penalty_iterator << ", best: ";
		//std::cerr << score_without_penalty << ", pen: " << actualPenalty << ", sum: ";
		//std::cerr << score_including_penalty << std::endl;
		
		if (  firstOptimization || score_including_penalty > best_score_including_penalty )
		{
			best_score_wo_penalty = score_without_penalty;
			best_score_including_penalty = score_including_penalty;
			std::copy( optimal_state_original.begin(), optimal_state_original.end(), optimal_state.begin()); 
		}
		theNaEManager->clear();
		firstOptimization = false;
	}
	gths.setAllStatesEnabled();

	
	//cerr<<"@@@ Optimal State: ";
	//for (i=0;i<numItems;i++)
	//{
	//	std::cerr<<(optimal_state[i]+1)<<" ";
	//}
	//std::cerr<< "score: " << best_score_wo_penalty << std::endl;
	if ( _outputNotice ) { std::cerr << " Optimal score following low resolution optimization: " << best_score_wo_penalty << std::endl;}

  double bestOScore = best_score_wo_penalty;											

	//cerr<<"hi5"<<endl;

	//---------- Orient Clique --------------------------------------//
	for (i=0; i<numItems; i++)
	{
    item[i]->initOrientation(*this);
        
		for (j=0; j<optimal_state[i];j++)
	   {
	   	item[i]->nextOrientation(*this); 
	   }
	}
	std::vector< int > optimalFlipState( numItems, -1 );
	for (i=0; i<numItems; i++)
	{
		bool badBumpBool;
		currScore[i] = item[i]->determineScore(*this,
	   	_dotBucket, _nBondCutoff, _probeRadius,
			_pmag, penalty, currBump[i],
			currHbond[i], badBumpBool);
		currBadBump[i] = badBumpBool;
		
		item[i]->setBestOrientation(
			item[i]->orientation(),
			currScore[i], currBadBump[i]);

		if (item[i]->canFlip())
		{
			//std::cerr << "about to track flip max score" << std::endl;
			item[i]->trackFlipMaxScore(item[i]->flipState(),
				best_score_wo_penalty, currBadBump[i]);
			optimalFlipState[ i ] = item[i]->flipState();
			if (item[i]->flipState() != 0)
			{
				item[i]->markFlipAtoms();
			}
			//std::cerr << "optimal flip state: " << i << " " << optimalFlipState[ i ] << std::endl;
		}
	}
	
	//determine the consequence for not flipping flipable-nodes
	for (int ii = 0; ii < numItems; ++ii)
	{
		if (! item[ii]->canFlip()) continue;
		//disable flip states opposite the best flip state.
		//std::cerr << "Determining best network state with sub-optimal flip for mover " << ii << std::endl;
		gths.setAllStatesEnabled();
		int const ii_num_states = gths.getNumStatesForNode( ii );
		item[ ii ]->initOrientation(*this );
		for (int jj = 0; jj < ii_num_states; ++jj)
		{
			if ( item[ ii ]->flipState() == optimalFlipState[ ii ] )			
			{
				//std::cerr << "Disabling State: " << ii << " " << jj << " " << item[ ii ]->flipState() << "==" << optimalFlipState[ ii ] << std::endl;
				gths.disableStateOnNode( ii, jj );	
			}
			item[ii]->nextOrientation(*this );
		}
		gths.setStateDisablingCompleteForNow();
		//std::cerr << "state disabling complete" << std::endl;
		if ( gths.anyNodeWithAllStatesDisabled() )
		{
			std::cerr << " Skipping iteration; node with all states disabled" << std::endl;
			continue;
		}
		
		NodeAndEdgeManager* theNaEManager = NodeAndEdgeManager::getInstance();
		theNaEManager->InitializeNetwork(gths);
		theNaEManager->computeOptimalNetworkConfiguration();
		//theNaEManager->bruteForceOriginalGraph();
		
		optimal_state_enabled = theNaEManager->getOptimalNetworkState();
		for (int jj = 0; jj < numItems; ++jj)
		{
			optimal_state_original[ jj ] = 
				gths.convertEnabledState2OriginalStateEnumerationOnNode(
				jj, optimal_state_enabled[ jj ] );
		}
		float score = theNaEManager->getScoreOfOptimalNetworkConfiguration();
		
		for (int jj = 0; jj < numItems; ++jj)
		{
			item[ jj ]->initOrientation(*this );
			for (int kk = 0; kk < optimal_state_original[ jj ]; ++kk)
			{
				item[ jj ]->nextOrientation(*this);
			}
		}
	
		//std::cerr << "Score According to gths: " << gths.getScoreForStateAssignment( optimal_state_original ) << std::endl;
		
		bool badBumpBool;
		currScore[ii] = item[ii]->determineScore(*this,
	   	_dotBucket, _nBondCutoff, _probeRadius,
			_pmag, penalty, currBump[ii],
			currHbond[ii], badBumpBool);
		currBadBump[ii] = badBumpBool;

		item[ii]->trackFlipMaxScore(item[ii]->flipState(), score, currBadBump[ii]);
		
		theNaEManager->clear();
	}
	gths.setAllStatesEnabled();

	//reset clique into optimal configuration
	for (i=0; i<numItems; i++)
	{
    item[i]->initOrientation(*this);
        
		for (j=0; j<optimal_state[i];j++)
	   {
	   	item[i]->nextOrientation(*this); 
	   }
	}
 
 // ---------------------------------------------------------------
  // optimization of any hires movers and get total score
  
  bool anyHighRes = false;
  int numInNonDefaultOrientations = 0;
  int cursor = 0;
  for (cursor = 0; cursor < numItems; cursor++) {
		if (item[cursor]->hasHires()) {
		anyHighRes = true;
      const int numHiResO = item[cursor]->numOrientations(Mover::HIGH_RES);
      // try each high-res orientation
      for (int ho = 0; ho < numHiResO; ho++) {
        if (item[cursor]->setOrientation(ho, *this, Mover::HIGH_RES)) {
          scoreThisO   = 0.0;
          worstPenalty = 0.0;
          if (_outputNotice && (! _showOrientScore) && _cliqueTicks){
            _os << ">";
          }
          bool tempbool;
          for (i = 0; i < numItems; i++) {
          	tempbool = currBadBump[i];
            currScore[i] = item[i]->determineScore(*this,
              _dotBucket, _nBondCutoff, _probeRadius,
              _pmag, penalty, currBump[i],
              currHbond[i], tempbool);
            currBadBump[i] = tempbool;
            if (penalty < worstPenalty) { worstPenalty = penalty; }
//            if (_outputNotice) {
//              if (_showOrientScore) {
//                _os << ">" <<item[i]->descr()
//                  << "[" << item[i]->orientation()+1;
//                if (i == cursor) { _os << "." << ho+1; }
//                _os << ":" << item[i]->describeOrientation()
//                  << "] bump=" << currBump[i]
//                  << ", HB=" << currHbond[i]
//                  << ((currBadBump[i]==TRUE)?", BADBUMP":"")
//                  << endl;
//              }
//              else if (_cliqueTicks) {
//                _os <<item[i]->orientation()+1;
//                if (i == cursor) { _os <<"."<<ho+1; }
//              	  _os <<" ";
//              }
//            }
            scoreThisO += currScore[i];
          }
          for (i = 0; i < numItems; i++) {
            if (item[i]->canFlip()) {
              item[i]->trackFlipMaxScore(item[i]->flipState(),
                scoreThisO, currBadBump[i]);
            }
          }
          
//          if ((ho == 0) || (scoreThisO+worstPenalty > bestOScore+prevWorstPenalty)) {
           if ((ho == 0) || (scoreThisO > bestOScore)) {

            bestOScore = scoreThisO;
//            prevWorstPenalty = worstPenalty;
            for (i = 0; i < numItems; i++) {
              item[i]->setBestOrientation(
                item[i]->orientation(Mover::HIGH_RES),
                currScore[i], currBadBump[i], Mover::HIGH_RES);
            }
          }
          if (_outputNotice) {
            if (_showOrientScore) {
              _os << "> adjustment score =" << scoreThisO
                << ", maxPenalty=" << worstPenalty
                << ", bestScore=" << bestOScore << endl;
            }
            else if (_cliqueTicks) {
              _os << " : " << bestOScore << "       \r";
            }
          }
        }
        else { break; } // could not set orientation (we are done)
      }
      // remember and position to the best
      const int bestHires = item[cursor]->bestOrientation(Mover::HIGH_RES);
      if (item[cursor]->orientation(Mover::HIGH_RES) != bestHires) {
        item[cursor]->setOrientation(bestHires, *this, Mover::HIGH_RES);
      }
      if (! item[cursor]->isDefaultO(bestHires, Mover::HIGH_RES)) {
        numInNonDefaultOrientations++;
      }
      bool tempbool;
      tempbool = currBadBump[cursor];
      currScore[cursor] = item[cursor]->determineScore(*this,
        _dotBucket, _nBondCutoff, _probeRadius,
        _pmag, penalty, currBump[cursor],
        currHbond[cursor], tempbool);
      currBadBump[cursor] = tempbool;
    }
    else {
      best = item[cursor]->bestOrientation();
      if (! item[cursor]->isDefaultO(best)) {
        numInNonDefaultOrientations++;
      }
    }
  }
  if (anyHighRes && _outputNotice ) { _os << std::endl;}
  // ---------------------------------------------------------------
  // tally the results
  if ( _outputNotice ) { _os << " Optimal score following high resolution, local optimization: " << bestOScore << endl;}
  
  return numInNonDefaultOrientations;
}


bool
AtomPositions::initializeCliqueMovers( 
	std::list< Mover* > const & clique, 
	std::vector< Mover* > & item,
	int const numItems
)
{
	std::copy( clique.begin(), clique.end(), item.begin() );
	for (int ii = 0; ii < numItems; ++ii)
	{
		if (! item[ii]->initOrientation(*this)) {
			cerr << endl << "ERROR: exhaustiveSearchOfClique: failed to initialize item - "
				<< item[ii]->descr() << endl;
			return false;
		}
		item[ii]->setBestOrientation(0, LowestMoverScore, FALSE);
	}
	return true;
}

void
AtomPositions::setNumStatesForNodes( 
	std::vector< Mover* > & item, 
	int const numItems,
	std::vector< int > & num_states,
	std::vector< std::vector< float > > & penalties
)
{
	double stateSpaceSize = 1;	
	//count number of states for movers
	for (int i=0; i<numItems; ++i)
	{
      item[i]->initOrientation(*this);
      int n_states=1;
      while (item[i]->nextOrientation(*this)) 
		{
         n_states++;
		}
		stateSpaceSize *= n_states;
	   num_states[ i ] = n_states;
	   penalties[ i ].resize( n_states );
	   std::fill( penalties[ i ].begin(), penalties[ i ].end(), float(0) );
   }
}

// ---------------------------------------------------------------
// insert user records about motion at the end of the header
void AtomPositions::describeChanges(std::list<PDBrec*>& records, 
									std::list<PDBrec*>::iterator& infoPtr, 
									std::vector<std::string>& notes) {
	if (notes.size() > 0) {
		PDBrec* separator = new PDBrec("USER  MOD -----------------------------------------------------------------");
		records.insert(infoPtr, separator);
		separator = new PDBrec("USER  MOD scores for adjustable sidechains, with \"set\" totals for H,N and Q");
		records.insert(infoPtr, separator);
		separator = new PDBrec("USER  MOD \"o\" means original, \"f\" means flipped, \"180deg\" is methyl default");
		records.insert(infoPtr, separator);
//		infoPtr.insertBefore(separator);
//		infoPtr.insertBefore(PDBrec(
//			"USER  MOD scores for adjustable sidechains, with \"set\" totals for H,N and Q"));
//		infoPtr.insertBefore(PDBrec(
//			"USER  MOD \"o\" means original, \"f\" means flipped, \"180deg\" is methyl default"));

		char fmtbuf[123];
		::sprintf(fmtbuf, "USER  MOD \"!\" flags a clash with an overlap of %.2fA or greater",
			_bad_bump_gap_cutoff);
		separator = new PDBrec(fmtbuf);
		records.insert(infoPtr, separator);

		::sprintf(fmtbuf, "USER  MOD flip categories: \"K\"=keep, \"C\"=clashes, \"X\"=uncertain, \"F\"=flip");
		separator = new PDBrec(fmtbuf);
		records.insert(infoPtr, separator);

		std::string s;

		std::vector<std::string>::iterator np = notes.begin();

		while(np != notes.end()) {
//			.next(s)) {	// insert notes in sorted order
//			PDBrec userRec(s.array());
			separator = new PDBrec((*np).c_str());
//                      cerr << separator->resname() << endl; 
			records.insert(infoPtr, separator);
			++np;
		}

		separator = new PDBrec("USER  MOD -----------------------------------------------------------------");
		records.insert(infoPtr, separator);
	}
}

// ---------------------------------------------------------------
void AtomPositions::manageMetals(const ResBlk& rblk) {
	std::multimap<std::string, PDBrec*> pdb = rblk.atomIt();
	std::string key, descr;

	std::multimap<std::string, PDBrec*>::const_iterator pdbit = pdb.begin();
	PDBrec* a = NULL;
	while(pdbit != pdb.end()) {
		key = pdbit->first;
		for (; pdbit != pdb.end() && pdbit->first == key; ++pdbit) {
			a = pdbit->second;
			if (a->hasProp(METALIC_ATOM)) {
				BumperPoint* bp = new BumperPoint(a->loc(), 0, 0, a->vdwRad());
				_excludePoints.push_front(bp);
			}
		}
	}
}

// ---------------------------------------------------------------
// create possible orientations for H atoms on waters (identified elsewhere)
// and store these Hs in the xyz table
void AtomPositions::generateWaterPhantomHs(std::list<PDBrec*>& waters) {
	char descrbuf[30];
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
					
					::sprintf(descrbuf, "%-3.3s%c%-3.3s%4d%c%-2.2s",
						(isAromRingAtom ? "/R/" : ""), rec->alt(),
						rec->resname(), rec->resno(),
						rec->insCode(), rec->chain());
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

float AtomPositions::determineScoreForMover(
  	Mover* mover,
  	std::vector< std::pair< AtomDescr, DotsForAtom * > > & atoms_to_score,
	double & penalty
)
{
	atoms_to_score_ptr_ = & atoms_to_score;
	scoreAtomsAndDotsInAtomsToScoreVector_ = true;
	bool tempbool;
	float bump, hbond;
		
	float score = mover->determineScore(*this,
		_dotBucket, _nBondCutoff, _probeRadius,
		_pmag, penalty, bump,
		hbond, tempbool);
   
	scoreAtomsAndDotsInAtomsToScoreVector_ = false;
	atoms_to_score_ptr_ = 0;
	return score;
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

