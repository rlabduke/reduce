// name: RotMethylSym.cpp
// author: J. Michael Word
// date written: 2/7/98
// purpose: Symmetry-aware implementation for RotMethyl

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

#include "../RotMethyl.h"
#include "../AtomPositions.h"

#define START_ANGLE 180.0
#define ROUGH_STEP   30
#define FINE_STEP     1

void RotMethyl::finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel,
						 AtomPositions &xyz, DotSphManager& dotBucket) {

	if (isComplete()) {

		// pre-build lists of bonded atoms

		const double approxNbondDistLimit = 3.0 + 0.5*nBondCutoff; // just a rule of thumb

		_rot.push_front(&_heavyAtom);
		for(std::list<PDBrec*>::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
			PDBrec* thisAtom = *alst;
			if (thisAtom->valid()) {
			  std::list<PDBrec*>* temp = new std::list<PDBrec*>();
			  std::list< std::pair<PDBrec*, Point3d> > neighborList = xyz.neighbors(thisAtom->loc(), thisAtom->covRad(),
				approxNbondDistLimit);
			  bondedList(thisAtom, neighborList, nBondCutoff, _rot, temp);
			  _bnded.push_back(temp);
            }
		}
		_rot.pop_front();
	}
}

int RotMethyl::makebumpers(NeighborList<BumperPoint*>& sym_bblks,
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
			sym_bblks.insert(bp);
			if (a->vdwRad() > maxVDWrad) { maxVDWrad = a->vdwRad(); }
		}
	}
	return an;
}
