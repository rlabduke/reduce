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

/*SJ what does f!=0 inside [] mean? f!=0 is a check that returns a 0 if f is equal to  0 and returns a 1 if f is any other number. So if f that is passed to these two functions is 0, then the first element of the array is changed/referenced, if f is any number other than a 0, the second element of the array is changed/referenced */
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

/// @brief Produce a list of atoms that are within a specified number of covalent-bond hops
/// from a specified atom.
///
/// Called by the finalize() method of Mover-derived classes to build a list of lists of atoms
/// bonded to the heavy atom or one of the hydrogens to be rotated.
///
/// @param [in] a The atom whose list of bonded atoms is to be returned.  This is
///				called for each atom in atmList by the finalize() method.
/// @param [in] nearby Atoms from the whole structure that are within a distance threshold that includes all
///				those that may be covalently bonded to the atom.
/// @param [in] nbnds How many bonds away before we stop counting (default is 3).
/// @param [in] atmList The heavy atom associated with a rotation along with all of the
///				hydrogens to be rotated.
/// @param [out] The function fills in this list of atoms from nearby that are bonded to "a".
void bondedList(const PDBrec& a, std::list< std::shared_ptr<PDBrec> >& nearby, int nbnds,
				std::list< std::shared_ptr<PDBrec> >& atmList, std::list< std::shared_ptr<PDBrec> >* bondedAtoms) {

	// Recursively tag all nearby atoms by how many covalent-bond hops they are away from the source.
	// First clear all the marks to 0.
	resetMarks(nearby);
	countBonds(a, nearby, 1, nbnds, atmList); // up to nbnds bonds away

	// Add all atoms that are from 1 to nbnds covalent hops from "a" to bondedAtoms
	std::shared_ptr<PDBrec> rec;
	for (std::list< std::shared_ptr<PDBrec> >::const_iterator it = nearby.begin(); it != nearby.end(); ++it) {
		rec = *it;
		if (rec->valid()) {
		  if (rec->mark() >= 1 && rec->mark() <= nbnds) {
			std::shared_ptr<PDBrec> temp = std::make_shared<PDBrec>(*rec);
			bondedAtoms->push_front(temp);
		  }
        }
	}
}

/// @brief Tag atoms by how many covalent bonds they are away from a source atom.
///
/// This function recursively tags nearby atoms by how far they are in covalent bond
/// hops from a source atom.  All atoms covalently bonded to the source are marked 1,
/// atoms covalently bonded to those are marked 2, and so forth.
/// This function has a side effect of tagging atoms in the nearby list with how many
/// covalent hops they are from the source atom.
/// @param [in] src Source atom to count bonds from.  For recursive calls, this
///				is called for all of the atoms that were marked at this level of recursion
///				so that we mark their neighbors.
/// @param [in] nearby Atoms from the whole structure that are close enough to be
///				potentially in the chains of covalently bonded atoms.  For recursive
///				calls, this is pared down to only atoms that have not yet been marked.
/// @param [in] distcount The current level of recursion.  The non-recursive caller should
///				set this to 1 to start counting direct neighbors.  Recursive calls
///				increase this by 1 every time to mark successive levels.
/// @param [in] maxcnt Terminate the recursion when we reach this many hops from the src.
/// @param [in] atmList The heavy atom associated with a rotation along with all of the
///				hydrogens to be rotated.
void countBonds(const PDBrec& src, const std::list< std::shared_ptr<PDBrec> >& nearby,
				int distcount, int maxcnt, std::list< std::shared_ptr<PDBrec> >& atmList) {

	// Keeps track of atoms that were marked at this level of recursion
	std::list< std::shared_ptr<PDBrec> > newlyMarked;
	// Keeps track of atoms that have not yet been marked so they should be passed to the next recursion level
	std::list< std::shared_ptr<PDBrec> > remainder;

	// Look at each of the atoms in the potentially-taggable list.  Split them into ones that are tagged
	// this iterations and ones to send on to future iterations.  Only consider records that have not been
	// tagged as invalid; invalid ones go on neither list.  Only check ones that are on the same alternate
	// conformations; others go onto neither list.
	std::shared_ptr<PDBrec> rec;
	for (std::list< std::shared_ptr<PDBrec> >::const_iterator targ = nearby.begin(); targ != nearby.end(); ++targ) {
		rec = *targ;
		if (rec->valid()) {
		  if ( (rec->mark() < 1 || rec->mark() > distcount)
			&& ! diffAltLoc(src, *rec) ) {

			// Atoms that are between -0.5 and 0.2 times the expected covalent-bond distances are
			// candidates to be bonded unless that bond is determined to be impossible.
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

	// So long as we have another valid hop count, re-run recursively on each atom that we just found
	// to be bonded with a one-larger count.  Only pass the set of atoms that have not yet been covalently
	// bonded at this or another earlier pass for consideration.
	if (distcount < maxcnt) {
		for(std::list< std::shared_ptr<PDBrec> >::const_iterator it = newlyMarked.begin(); it != newlyMarked.end(); ++it) {
			countBonds(**it, remainder, distcount+1, maxcnt, atmList);
		}
	}
}

// initialize the markers we use to determine bonding patterns
void resetMarks(std::list< std::shared_ptr<PDBrec> >& lst) {
	for (std::list< std::shared_ptr<PDBrec> >::const_iterator it = lst.begin(); it != lst.end(); ++it) {
		(*it)->mark(0);
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

int withinCovalentDist(const PDBrec& p, const PDBrec& q, double offset) {
   double lim = p.covRad() + q.covRad() + offset;

   return distanceSquared(p.loc(), q.loc()) <= (lim*lim);
}

/// @brief Determine whether it is impossible for a covalent bond to form between two atoms.
///
/// This is decided in the context of a list of one or more hydrogen atoms to be attached to a
/// heavy atom in a residue.  One or both of the two atoms may be in the list of
/// hydrogen atoms attached to a heavy atom in a residue.  The atoms may be in the same or
/// in different residues.
/// @param [in] src First atom
/// @param [in] targ Second atom
/// @param [in] atmList The heavy atom associated with a rotation along with all of the
///				hydrogens to be rotated.  Exactly one of these atoms is not hydrogen.
bool impossibleCovalent(const PDBrec& src, const PDBrec& targ, std::list< std::shared_ptr<PDBrec> >& atmList)
{
   // If both of the atoms are hydrogens, then they cannot be covalently bonded to one another.
   if (src.isHydrogen() && targ.isHydrogen()) { return TRUE; }

   // Handle the case where exactly one of the atoms is hydrogen.
   else if (src.isHydrogen() || targ.isHydrogen()) {
	  // See whether each of the atoms is in the list that includes the heavy atom along with
	  // the hydrogens to be bonded to it.
      bool srcInList = FALSE, targInList = FALSE;
      for(std::list< std::shared_ptr<PDBrec> >::const_iterator it = atmList.begin(); it != atmList.end(); ++it) {
	      if (**it == src)  {srcInList  = TRUE; }
	 else if (**it == targ) {targInList = TRUE; }
	  }

	  // If both of them are in this list (and we already know that exactly one of them is
	  // a hydrogen) then it is not impossible for them to be bonded.  However, if one of
	  // them is a hydrogen in the list and the other is not in the list, then it is not
	  // possible for them to be covalently bonded because we already know that all of the
	  // bonds have been satisfied by this rotatable group.
      if (srcInList && targInList) { return FALSE; }
      else if (srcInList && src.isHydrogen())   { return TRUE; }
      else if (targInList && targ.isHydrogen()) { return TRUE; }
   }

   // If the two atoms are not in the same residue and at least one of them is a
   // hydrogen, then they cannot be covalently bonded.
   return (! sameres(src, targ)) && (src.isHydrogen() || targ.isHydrogen());
}

bool foundInList(const PDBrec& a, const std::list< std::shared_ptr<PDBrec> >& lst) {
	for (std::list< std::shared_ptr<PDBrec> >::const_iterator it = lst.begin(); it != lst.end(); ++it) {
		if (a == **it) { return TRUE; }
	}
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

bool annularDots(const Point3d& dot, const PDBrec& src, const PDBrec& targ, float probeRadius) {
   return dot2srcCenter(dot, src, targ) >
      kissEdge2bullsEye(src.vdwRad(), targ.vdwRad(), probeRadius);
}

double dot2srcCenter(const Point3d& dot, const PDBrec& src, const PDBrec& targ) {

   const Point3d src2targVec = (targ.loc() - src.loc()).scaleTo(src.vdwRad());
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
