// name: Mover.C
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

#include "Mover.h"

// ResBlk.h is included to define sameres()
#include "ResBlk.h"

void Mover::resetFlipMaxScore(int f) {
     _flipMaxScore[f!=0] = LowestMoverScore;
   _flipMaxBadBump[f!=0] = FALSE;
}

void Mover::trackFlipMaxScore(int f, double val, bool hasBadBump) {
   if ( _flipMaxScore[f!=0] < val) {
        _flipMaxScore[f!=0] = val;
      _flipMaxBadBump[f!=0] = hasBadBump;
   }
}

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

bool visableAltConf(const PDBrec& a, bool onlyA) {
   const char aalt = a.alt();
   bool returnvalue = (aalt == ' ' || aalt == 'A' || aalt == 'a' || !onlyA);
	//std::cerr << "visiableAltConf: " << a.recName() << " "  << returnvalue << std::endl;
	return returnvalue;
}

bool interactingConfs(const PDBrec& a, const PDBrec& b, bool onlyA) {
   return ( (!visableAltConf(a, onlyA)) || (!visableAltConf(b, onlyA)) )
      ? FALSE
      : (! diffAltLoc(a, b));
}

// are a and b alternative locations for the same residue?
bool diffAltLoc(const PDBrec& a, const PDBrec& b) {
   //return a.resno() == b.resno() && a.insCode() == b.insCode()
   return a.insCode() == b.insCode()
       && strcmp(a.chain(), b.chain()) == 0
       && a.alt() != b.alt() && a.alt() != ' ' && b.alt() != ' ';
}

int withinCovalentDist(std::pair<PDBrec*, Point3d> p, std::pair<PDBrec*, Point3d> q, double offset) {
   double lim = p.first->covRad() + q.first->covRad() + offset;

   return distanceSquared(p.second, q.second) <= (lim*lim);
}

bool impossibleCovalent(const PDBrec& src, const PDBrec& targ, std::list<PDBrec*>& atmList) {
   if (src.isHydrogen() && targ.isHydrogen()) { return TRUE; }

   else if (src.isHydrogen() || targ.isHydrogen()) {
      bool srcInList = FALSE, targInList = FALSE;
      for(std::list<PDBrec*>::const_iterator it = atmList.begin(); it != atmList.end(); ++it) {
	      if (**it == src)  {srcInList  = TRUE; }
	 else if (**it == targ) {targInList = TRUE; }
      }

      if (srcInList && targInList) { return FALSE; }
      else if (srcInList && src.isHydrogen())   { return TRUE; }
      else if (targInList && targ.isHydrogen()) { return TRUE; }
   }

   return (! sameres(src, targ)) && (src.isHydrogen() || targ.isHydrogen());
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

// what is the gap between the VDW radii of atom a at point aloc
// and atom b at point bloc?

double vdwGap(const PDBrec& a, const Point3d& aloc,
              const PDBrec& b, const Point3d& bloc) {
   double lim = a.vdwRad() + b.vdwRad();

   return distance2(aloc, bloc) - lim;
}

// functions used to restrict annular rings of good dots around clashes

bool annularDots(const Point3d& dot, const PDBrec& src, const PDBrec& targ, Point3d targ_loc, float probeRadius) {
   return dot2srcCenter(dot, src, targ_loc) >
      kissEdge2bullsEye(src.vdwRad(), targ.vdwRad(), probeRadius);
}

double dot2srcCenter(const Point3d& dot, const PDBrec& src, Point3d targ_loc) {

   const Point3d src2targVec = (targ_loc - src.loc()).scaleTo(src.vdwRad());
   const Point3d srcSurfacePoint = src2targVec + src.loc();
   return distance2(dot, srcSurfacePoint);
}

double kissEdge2bullsEye(float ra, float rb, float rp) {
   return 2.0*ra*sqrt(rb*rp/((ra+rb)*(ra+rp)));
}

// force sorting of Seq<MoverPtr> by descr()
// by duplicating some merge-sort code

void mpSeqSort(std::list<MoverPtr>& x) {
	x.reverse();
}

/*
std::list<MoverPtr> mpSeqMerge(const std::list<MoverPtr>& x_list, const std::list<MoverPtr>& y_list) {
	std::list<MoverPtr> r;
	std::list<MoverPtr>::const_iterator x = x_list.begin();
	std::list<MoverPtr>::const_iterator y = y_list.begin();

   while ((x != x_list.end()) && (y != y_list.end())) {
      if ( (*x) && (*y) && ((*x)->descr() < (*y)->descr())) {
         r.push_front(*x);
         ++x;
      } else {
         r.push_front(*y);
         ++y;
      }
   }
   while (x != x_list.end()) {
      r.push_front(*x);
      ++x;
   }
   while (y != y_list.end()) {
      r.push_front(*y);
      ++y;
   }
   r.reverse();
   return r;
}

void mpSeqSplit(std::list<MoverPtr> x_list, std::list<MoverPtr>& y, std::list<MoverPtr>& z) {
	std::list<MoverPtr>::const_iterator x = x_list.begin();
	while(x != x_list.end()) {
		y.push_front(*x);
		++x;
		if (x != x_list.end()) {
			z.push_front(*x);
			++x;
		}
	}
}

std::list<MoverPtr> mpSeqSort(const std::list<MoverPtr>& x) {
   if (x.empty() || (x.size() == 1)) { return x; }
   std::list<MoverPtr> p, q;
   mpSeqSplit(x, p, q);
   return mpSeqMerge(mpSeqSort(p), mpSeqSort(q));
}

bool mpLess(const MoverPtr& x, const MoverPtr& y) {
	return x->descr() < y->descr();
}
*/
