// name: MoverSym.cpp
// author: J. Michael Word
// date written: 2/7/98
// purpose: Symmetry-aware utility routines for bonding and geometric relationships

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

void bondedList(PDBrec *a, std::list< std::pair<PDBrec*, Point3d> >& nearby, int nbnds,
				std::list<PDBrec*>& atmList, std::list<PDBrec*>* bondedAtoms) {
	resetMarks(nearby);
	countBonds(std::make_pair(a, a->loc()), nearby, 1, nbnds, atmList); // up to nbnds bonds away

	PDBrec* rec = NULL;
	for (std::list< std::pair<PDBrec*, Point3d> >::const_iterator it = nearby.begin(); it != nearby.end(); ++it) {
		rec = it->first;
		if (rec->valid()) {
		  if (rec->mark() >= 1 && rec->mark() <= nbnds) {
			if (distanceSquared(rec->loc(),it->second) < 0.000001){ //Vishal: use TOL
				PDBrec* temp = new PDBrec(*rec);
				bondedAtoms->push_front(temp);
			}
		  }
        }
	}
}

void countBonds(std::pair<PDBrec*,Point3d> src, const std::list< std::pair<PDBrec*, Point3d> >& nearby,
				int distcount, int maxcnt, std::list<PDBrec*>& atmList) {

	std::list< std::pair<PDBrec*, Point3d> > newlyMarked;
	std::list< std::pair<PDBrec*, Point3d> > remainder;

	PDBrec* rec = NULL;
	for (std::list< std::pair<PDBrec*, Point3d> >::const_iterator targ = nearby.begin(); targ != nearby.end(); ++targ) {
		rec = targ->first;
		bool is_sym= (distanceSquared(rec->loc(),targ->second) > 0.00001); //Vishal: use TOL
		if (rec->valid()) {
		  if ( (rec->mark() < 1 || rec->mark() > distcount)
			&& ! diffAltLoc(*(src.first), *rec) ) {

			bool isnear   = withinCovalentDist(src, *targ,  0.2);
			bool tooclose = withinCovalentDist(src, *targ, -0.5);

			if (isnear && ! (is_sym || tooclose || impossibleCovalent(*(src.first), *rec, atmList))) {
				rec->mark(distcount);
				newlyMarked.push_front(*targ);
			}
			else {
				remainder.push_front(*targ);
			}
		  }
		}
	}
	if (distcount < maxcnt) {
		for(std::list< std::pair<PDBrec*, Point3d> >::const_iterator it = newlyMarked.begin(); it != newlyMarked.end(); ++it) {
			countBonds(*it, remainder, distcount+1, maxcnt, atmList);
		}
	}
}

// initialize the markers we use to determine bonding patterns
void resetMarks(std::list< std::pair<PDBrec*,Point3d> >& lst) {
	for (std::list< std::pair<PDBrec*, Point3d> >::const_iterator it = lst.begin(); it != lst.end(); ++it) {
		it->first->mark(0);
	}
}

int withinCovalentDist(std::pair<PDBrec*, Point3d> p, std::pair<PDBrec*, Point3d> q, double offset) {
   double lim = p.first->covRad() + q.first->covRad() + offset;

   return distanceSquared(p.second, q.second) <= (lim*lim);
}

bool foundInList(const PDBrec& a, Point3d p, const std::list<PDBrec*>& lst) {
	//PDBrec a_copy=a;
	//a_copy.newPDBrecRep();
	//a_copy.loc(p);
	if (distanceSquared(a.loc(),p) > 0.000001)
		return FALSE;
	for (std::list<PDBrec*>::const_iterator it = lst.begin(); it != lst.end(); ++it) {
		if (a == **it){ return TRUE; }
	}
	//a_copy.delPDBrecRep();
	return FALSE;
}

bool annularDots(const Point3d& dot, const PDBrec& src, const PDBrec& targ, Point3d targ_loc, float probeRadius) {
   return dot2srcCenter(dot, src, targ_loc) >
      kissEdge2bullsEye(src.vdwRad(), targ.vdwRad(), probeRadius);
}

double dot2srcCenter(const Point3d& dot, const PDBrec& src, Point3d targ_loc) {

   const Point3d src2targVec = (targ_loc - src.loc()).scaleTo(src.vdwRad());
   const Point3d srcSurfacePoint = src2targVec + src.loc();
   return distance2(dot, srcSurfacePoint);
}
