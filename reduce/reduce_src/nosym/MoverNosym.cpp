// name: MoverNosym.cpp
// author: J. Michael Word
// date written: 2/7/98
// purpose: utility routines for bonding and geometric relationships

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#include "../Mover.h"

// ResBlk.h is included to define sameres()
#include "../ResBlk.h"

void bondedList(const PDBrec& a, std::list<PDBrec*>& nearby, int nbnds,
				std::list<PDBrec*>& atmList, std::list<PDBrec*>* bondedAtoms) {
	resetMarks(nearby);
	countBonds(a, nearby, 1, nbnds, atmList); // up to nbnds bonds away

	PDBrec* rec = NULL;
	for (std::list<PDBrec*>::const_iterator it = nearby.begin(); it != nearby.end(); ++it) {
		rec = *it;
		if (rec->valid()) {
		  if (rec->mark() >= 1 && rec->mark() <= nbnds) {
			PDBrec* temp = new PDBrec(*rec);
			bondedAtoms->push_front(temp);
		  }
        }
	}
}

void countBonds(const PDBrec& src, const std::list<PDBrec*>& nearby,
				int distcount, int maxcnt, std::list<PDBrec*>& atmList) {

	std::list<PDBrec*> newlyMarked;
	std::list<PDBrec*> remainder;

	PDBrec* rec = NULL;
	for (std::list<PDBrec*>::const_iterator targ = nearby.begin(); targ != nearby.end(); ++targ) {
		rec = *targ;
		if (rec->valid()) {
		  if ( (rec->mark() < 1 || rec->mark() > distcount)
			&& ! diffAltLoc(src, *rec) ) {

			bool isnear   = withinCovalentDist(src, *rec,  0.2);
			bool tooclose = withinCovalentDist(src, *rec, -0.5);

			if (isnear && ! (tooclose || impossibleCovalent(src, *rec, atmList))) {
				rec->mark(distcount);
				newlyMarked.push_front(rec);
			}
			else {
				remainder.push_front(rec);
			}
		  }
		}
	}
	if (distcount < maxcnt) {
		for(std::list<PDBrec*>::const_iterator it = newlyMarked.begin(); it != newlyMarked.end(); ++it) {
			countBonds(**it, remainder, distcount+1, maxcnt, atmList);
		}
	}
}

// initialize the markers we use to determine bonding patterns
void resetMarks(std::list<PDBrec*>& lst) {
	for (std::list<PDBrec*>::const_iterator it = lst.begin(); it != lst.end(); ++it) {
		(*it)->mark(0);
	}
}

int withinCovalentDist(const PDBrec& p, const PDBrec& q, double offset) {
   double lim = p.covRad() + q.covRad() + offset;

   return distanceSquared(p.loc(), q.loc()) <= (lim*lim);
}

bool foundInList(const PDBrec& a, const std::list<PDBrec*>& lst) {
	for (std::list<PDBrec*>::const_iterator it = lst.begin(); it != lst.end(); ++it) {
		if (a == **it) { return TRUE; }
	}
	return FALSE;
}

bool annularDots(const Point3d& dot, const PDBrec& src, const PDBrec& targ, float probeRadius) {
   return dot2srcCenter(dot, src, targ) >
      kissEdge2bullsEye(src.vdwRad(), targ.vdwRad(), probeRadius);
}

double dot2srcCenter(const Point3d& dot, const PDBrec& src, const PDBrec& targ) {

   const Point3d src2targVec = (targ.loc() - src.loc()).scaleTo(src.vdwRad());
   const Point3d srcSurfacePoint = src2targVec + src.loc();
   return distance2(dot, srcSurfacePoint);
}
