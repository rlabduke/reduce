// name: RotAromMethylNosym.cpp
// author: J. Michael Word, modified by Aram Han
// date written: 2/7/98, modified 8/13/12
// purpose: Implementation for RotAromMethyl

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
#include <sstream>
using std::cout;
#endif

#include "../RotAromMethyl.h"
#include "../AtomPositions.h"

void RotAromMethyl::finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel,
						 AtomPositions &xyz, DotSphManager& dotBucket) {

	if (isComplete()) {

		// pre-build lists of bonded atoms

		const double approxNbondDistLimit = 3.0 + 0.5*nBondCutoff; // just a rule of thumb

		_rot.push_front(&_heavyAtom);
		for(std::list<PDBrec*>::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
			const PDBrec* thisAtom = *alst;
			std::list<PDBrec*>* temp = new std::list<PDBrec*>();
			std::list<PDBrec*> neighborList = xyz.neighbors(thisAtom->loc(), thisAtom->covRad(),
				approxNbondDistLimit);
			bondedList(*thisAtom, neighborList, nBondCutoff, _rot, temp);
			_bnded.push_back(temp);
		}
		_rot.pop_front();
	}
}

int RotAromMethyl::makebumpers(std::multimap<LocBlk, BumperPoint*>& bblks,
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
