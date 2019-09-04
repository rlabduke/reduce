// name: FlipMemo.C
// author: J. Michael Word
// date written: 2/7/98
// purpose: Implementation for FlipMemo

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
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#else
#include <cctype>
#include <cmath>
#include <cstdlib>
using std::exit;
using std::toupper;
#endif

#include "FlipMemo.h"
#include "AtomConn.h"
#include "AtomPositions.h"

extern bool UseNuclearDistances; //defined in reduce.cpp JJH
extern bool GenerateFinalFlip; // SJ - defined in reduce.cpp for checking if flips are being generated for scoring or final PDB file

const char* FlipMemo::_pointName[FMnumFlipRes+1]
                                [FMmaxAtomSlots+1] = {
{"HIS", " ND1", " CD2", " CE1", " NE2",         // all atoms needed to det.
        " HD1", " HD2", " HE1", " HE2", " CG",  // locs together at front
	" CB", " CA", " HA", " N", " C", "1HB", "2HB", // PDBv2.3 - USEOLD
	" HD1", " HD2", " HE1", " HE2"},
{"HIS", " ND1", " CD2", " CE1", " NE2",         // PDBv3.0 - USENEW
        " HD1", " HD2", " HE1", " HE2", " CG",
        " CB", " CA", " HA", " N", " C", " HB2", " HB3",
        " HD1", " HD2", " HE1", " HE2"},
{"ASN", " OD1", " ND2", "1HD2", "2HD2", " CG",        // Not Xplor = PDBv2.3 - USEOLD
        " CB", " CA", " HA", " N", " C", "1HB", "2HB",
	"1HD2", "2HD2"},
{"ASN", " OD1", " ND2", "HD21", "HD22", " CG",        // Xplor
        " CB", " CA", " HA", " N", " C", "1HB", "2HB",
	"HD21", "HD22"},
{"ASN", " OD1", " ND2", "HD21", "HD22", " CG",        // PDBv3.0 - USENEW
        " CB", " CA", " HA", " N", " C", " HB2", " HB3",
        "HD21", "HD22"},
{"GLN", " OE1", " NE2", "1HE2", "2HE2", " CD",        // Not Xplor = PDBv2.3 - USEOLD
        " CG", " CB", " CA", "1HB", "2HB", "1HG", "2HG",
	"1HE2", "2HE2"},
{"GLN", " OE1", " NE2", "HE21", "HE22", " CD",        // Xplor
        " CG", " CB", " CA", "1HB", "2HB", "1HG", "2HG",
	"HE21", "HE22"},
{"GLN", " OE1", " NE2", "HE21", "HE22", " CD",        // PDBv3.0 - USENEW
        " CG", " CB", " CA", " HB2", " HB3", " HG2", " HG3",
        "HE21", "HE22"},
{NULL, NULL}
};

// SJ - 09/15/2015 - Index of atoms in the above _pointName array (and hence the _origLoc array and the _wrkAtom array) involved in the three step flip. First number is the number of atoms
const int FlipMemo::_dockAtomIndex[FMnumFlipRes][FMmaxAtomSlots+1] = {
    {13, 10/*CB*/,9/*CG*/,3/*CE1*/,4/*NE2*/,11/*CA*/,1/*ND1*/,2/*CD2*/,5/*HD1*/,6/*HD2*/,7/*HE1*/,8/*HE2*/,15/*1HB*/,16/*2HB*/}, // HIS - PDBv2.3 - USEOLD
    {13, 10/*CB*/,9/*CG*/,3/*CE1*/,4/*NE2*/,11/*CA*/,1/*ND1*/,2/*CD2*/,5/*HD1*/,6/*HD2*/,7/*HE1*/,8/*HE2*/,15/*HB2*/,16/*HB3*/}, // HIS - PDBv3.0 - USENEW
    {9, 6/*CB*/,5/*CG*/,1/*OD1*/,2/*ND2*/,7/*CA*/,3/*1HD2*/,4/*2HD2*/,11/*1HB*/,12/*2HB*/}, // ASN - Not Xplor = PDBv2.3 - USEOLD
    {9, 6/*CB*/,5/*CG*/,1/*OD1*/,2/*ND2*/,7/*CA*/,3/*HD21*/,4/*HD22*/,11/*1HB*/,12/*2HB*/}, // ASN - Xplor
    {9, 6/*CB*/,5/*CG*/,1/*OD1*/,2/*ND2*/,7/*CA*/,3/*HD21*/,4/*HD22*/,11/*HB2*/,12/*HB3*/}, // ASN - PDBv3.0 - USENEW
    {12, 6/*CG*/,5/*CD*/,1/*OE1*/,2/*NE2*/,8/*CA*/,3/*1HE2*/,4/*2HE2*/,7/*CB*/,9/*1HB*/,10/*2HB*/,11/*1HG*/,12/*2HG*/}, //GLN - Not Xplor = PDBv2.3 - USEOLD
    {12, 6/*CG*/,5/*CD*/,1/*OE1*/,2/*NE2*/,8/*CA*/,3/*HE21*/,4/*HE22*/,7/*CB*/,9/*1HB*/,10/*2HB*/,11/*1HG*/,12/*2HG*/}, //GLN - Xplor
    {12, 6/*CG*/,5/*CD*/,1/*OE1*/,2/*NE2*/,8/*CA*/,3/*HE21*/,4/*HE22*/,7/*CB*/,9/*HB2*/,10/*HB3*/,11/*HG2*/,12/*HG3*/} //GLN - PDBv3.0 - USENEW
};

const ResFlipInfo_t FlipMemo::_resFlip[FMnumFlipRes+1] = {
// --- rname, fromO,numO,numXO, fromPP,numPP, from3B,numScore,numBmpr,numHeavy,
// ---        fromScat,numScat, numPnts, flags

#ifdef AROMATICS_ACCEPT_HBONDS
// {"HIS", 0,6,8, 0,4, 0,9,9,4, 1,9, 20,           0},
 {"HIS", 0,6,8, 0,4, 0,9,9,4, 1,9, 20, USEOLDNAMES}, // SJ - 09/16/2015 changed the last argument from 0 to USEOLDNAMES - was a bug - TODO have to check if it messes anything else up
 {"HIS", 0,6,8, 0,4, 0,9,9,4, 1,9, 20, USENEWNAMES},
#else
 //{"HIS", 0,6,8, 0,4, 0,9,9,4, 1,9, 20,           0},
 {"HIS", 0,6,8, 0,4, 0,8,8,4, 1,9, 20, USEOLDNAMES}, // SJ - 09/16/2015 changed the last argument from 0 to USEOLDNAMES - was a bug - TODO have to check if it messes anything else up
 {"HIS", 0,6,8, 0,4, 0,9,9,4, 1,9, 20, USENEWNAMES},
#endif
 {"ASN", 8,2,2, 4,2, 9,5,4,2, 1,5, 14, USEOLDNAMES},
 {"ASN", 8,2,2, 4,2, 9,5,4,2, 1,5, 14,   XPLORNAME},
 {"ASN", 8,2,2, 4,2, 9,5,4,2, 1,5, 14, USENEWNAMES},
 {"GLN", 8,2,2, 6,2,14,5,4,2, 1,5, 14, USEOLDNAMES},
 {"GLN", 8,2,2, 6,2,14,5,4,2, 1,5, 14,   XPLORNAME},
 {"GLN", 8,2,2, 6,2,14,5,4,2, 1,5, 14, USENEWNAMES},
 { NULL, 0,0,0, 0,0, 0,0,0,0, 0,0,  0,           0}
};

const ProtonLoc_t FlipMemo::_protLoc[FMnumNewLocs+1] = {
 {"HIS", " HD1", 4, 17, 2, 9, 4, 0.86, 1.02,   0.0,   0.0},
 {"HIS", " HD2", 4, 18, 1, 3, 9, 0.93, 1.08,   0.0,   0.0},
 {"HIS", " HE1", 4, 19, 4, 2, 3, 0.93, 1.08,   0.0,   0.0},
 {"HIS", " HE2", 4, 20, 3, 4, 1, 0.86, 1.02,   0.0,   0.0},
 {"ASN", "1HD2", 3, 13, 1, 5, 2, 0.86, 1.02, 120.0, 180.0},
 {"ASN", "2HD2", 3, 14, 1, 5, 2, 0.86, 1.02, 120.0,   0.0},
 {"GLN", "1HE2", 3, 13, 1, 5, 2, 0.86, 1.02, 120.0, 180.0},
 {"GLN", "2HE2", 3, 14, 1, 5, 2, 0.86, 1.02, 120.0,   0.0},
 { NULL,   NULL, 0,  0, 0, 0, 0, 0.0,  0.0,    0.0,   0.0}
};

const BondLimits_t FlipMemo::_bondedLimits[FnumScoreableAtoms] = {
 { 1,  3, 7,12 }, // HIS...
 { 2,  3, 7,12 },
 { 3,  3, 7, 9 },
 { 4,  3, 7, 9 },
 { 5,  1, 3, 7 },
 { 6,  1, 3, 7 },
 { 7,  1, 3, 7 },
 { 8,  1, 3, 7 },
 { 9,  3,10,15 },
 { 1,  1, 3, 8 }, // ASN
 { 2,  3, 5, 8 },
 { 3,  1, 3, 5 },
 { 4,  1, 3, 5 },
 { 5,  3, 8,11 },
 { 1,  1, 3, 8 }, // GLN
 { 2,  3, 5, 8 },
 { 3,  1, 3, 5 },
 { 4,  1, 3, 5 },
 { 5,  3, 8,11 },
};
const int FlipMemo::_bondedSet[FnumScoreableAtoms]
                              [FMmaxBondedAtoms] = {
 { 3, 5, 9, 2, 4, 7,10, 6, 8,11,15,16 }, // HIS...
 { 4, 6, 9, 1, 3, 8,10, 5, 7,11,15,16 },
 { 1, 4, 7, 2, 5, 8, 9, 6,10 },
 { 2, 3, 8, 1, 6, 7, 9, 5,10 },
 { 1, 3, 9, 2, 4, 7,10 },
 { 2, 4, 9, 1, 3, 8,10 },
 { 3, 1, 4, 2, 5, 8, 9 },
 { 4, 2, 3, 1, 6, 7, 9 },
 { 1, 2,10, 3, 4, 5, 6,11,15,16, 7, 8, 12, 13, 14 },
 { 5, 2, 6, 3, 4, 7,11,12 },             // ASN...
 { 3, 4, 5, 1, 6, 7,11,12 },
 { 2, 4, 5, 1, 6 },
 { 2, 3, 5, 1, 6 },
 { 1, 2, 6, 3, 4, 7,11,12, 8, 9,10 },
 { 5, 2, 6, 3, 4, 7,11,12 },             // GLN...
 { 3, 4, 5, 1, 6, 7,11,12 },
 { 2, 4, 5, 1, 6 },
 { 2, 3, 5, 1, 6 },
 { 1, 2, 6, 3, 4, 7,11,12, 8, 9,10 },
};

const int FlipMemo::_atomOrient[FMnumFlipOrient]
                               [FMmaxAtomSlots+1] = { // SJ this is how the coordinates of the atoms are exchanged for flips
#ifdef AROMATICS_ACCEPT_HBONDS
 {0, 1, 2, 3, 4, 0, 6, 7, 8, 9}, // HIS... // SJ - when the number is 0 (apart from the first number) it means that atom does not exist - this is for HIS only with different protonation states
 {1, 1, 2, 3, 4, 5, 6, 7, 0, 9},
 {2, 1, 2, 3, 4, 5, 6, 7, 8, 9},
 {3, 2, 1, 4, 3, 0,18,19,20, 9}, // SJ - flipped
 {4, 2, 1, 4, 3,17,18,19, 0, 9}, // flipped
 {5, 2, 1, 4, 3,17,18,19,20, 9}, // flipped
 {6, 1, 2, 3, 4, 0, 6, 7, 0, 9}, // rarely used double deprotonated states
 {7, 2, 1, 4, 3, 0,18,19, 0, 9}, // flipped for rarely used double deprotonated states
#else
 {0, 1, 2, 3, 4, 0, 6, 7, 8}, // HIS...
 {1, 1, 2, 3, 4, 5, 6, 7, 0},
 {2, 1, 2, 3, 4, 5, 6, 7, 8},
 {3, 2, 1, 4, 3, 0,18,19,20}, // flipped
 {4, 2, 1, 4, 3,17,18,19, 0}, // flipped
 {5, 2, 1, 4, 3,17,18,19,20}, // flipped
 {6, 1, 2, 3, 4, 0, 6, 7, 0}, // rarely used double deprotonated states
 {7, 2, 1, 4, 3, 0,18,19, 0}, // flipped for rarely used double deprotonated states
#endif
 {0, 1, 2, 3, 4, 5},          // ASN, GLN... // SJ - unflipped
 {1, 2, 1,13,14, 5}                          // SJ - flipped - don't know what the 5 means at the end.
};

const int FlipMemo::_atomDAflags[FMnumFlipOrient]
                                [FMmaxAtomSlots+1] = {
#ifdef AROMATICS_ACCEPT_HBONDS
 {0,-1,-1,-1,-1, 0, 0, 0,+1,-1}, // HIS...
 {1,-1,-1,-1,-1,+1, 0, 0, 0,-1},
 {2, 0, 0, 0, 0,+1, 0, 0,+1, 0}, // set orientation uses this
 {3,-1,-1,-1,-1, 0, 0, 0,+1,-1}, // only to set nitrogen and carbon D/A status
 {4,-1,-1,-1,-1,+1, 0, 0, 0,-1},
 {5, 0, 0, 0, 0,+1, 0, 0,+1, 0},
 {6,-1, 0, 0,-1, 0, 0, 0, 0, 0}, // rarely used double deprotonated states
 {7,-1, 0, 0,-1, 0, 0, 0, 0, 0},
#else
 {0,-1, 0, 0, 0, 0, 0, 0,+1, 0}, // HIS...
 {1, 0, 0, 0,-1,+1, 0, 0, 0, 0},
 {2, 0, 0, 0, 0,+1, 0, 0,+1, 0}, // set orientation uses this
 {3,-1, 0, 0, 0, 0, 0, 0,+1, 0}, // only to set nitrogen D/A status
 {4, 0, 0, 0,-1,+1, 0, 0, 0, 0},
 {5, 0, 0, 0, 0,+1, 0, 0,+1, 0},
 {6,-1, 0, 0,-1, 0, 0, 0, 0, 0}, // rarely used double deprotonated states
 {7,-1, 0, 0,-1, 0, 0, 0, 0, 0},
#endif
 {0,-1, 0,+1,+1, 0},             // ASN, GLN...
 {1,-1, 0,+1,+1, 0}
};

const char* FlipMemo::_orientDescr[FMnumFlipOrient] = {
 "     no HD1", // HIS...
 "     no HE2",
 "    +bothHN",
 "FLIP no HD1",
 "FLIP no HE2",
 "FLIP+bothHN",
 "    - no HN", // rarely used double deprotonated states
 "FLIP- no HN",
 "      amide", // ASN, GLN...
 "FLIP  amide"

#ifdef BADXXX
 "       w/o  HD1", // HIS...
 "       w/o  HE2",
 "     + both HN ",
 "FLIP   w/o  HD1",
 "FLIP   w/o  HE2",
 "FLIP + both HN ",
 "     -   no HN ", // rarely used double deprotonated states
 "FLIP -   no HN ",
 "       sc amide", // ASN, GLN...
 "FLIP   sc amide"
#endif
};

const double FlipMemo::_orientationPenalty[FMnumFlipOrient] = {
 0.00, // HIS...
 0.00,
 0.05,
 0.50,
 0.50,
 0.55,
 1.00, // double deprotonated HIS - tried when extendOrientations(TRUE)
 1.50,
 0.00, // ASN, GLN...
 0.50
};

const int FlipMemo::_oFlipState[FMnumFlipOrient] = {
 0, // HIS...
 0,
 0,
 1,
 1,
 1,
 0, // rarely used double deprotonated states
 1,
 0, // ASN, GLN...
 1
};

FlipMemo::FlipMemo(const char *resname, bool useXplorNames, bool useOldNames, bool bbModel) {
   _isComplete = FALSE;
   _extendO    = FALSE;
   _fromO      = 0;
   _numO       = 0;

   // look for residue type by comparing residue names
   // and comparing the flags (use XPLOR names?)
   // if not found, the one we stop on is blank

    for(_resType = 0; _resFlip[_resType].rname; _resType++) {
      if ( (strcmp(_resFlip[_resType].rname, resname) == 0)
      && !(((_resFlip[_resType].flags & USEOLDNAMES) && ! useOldNames)
        || ((_resFlip[_resType].flags & XPLORNAME)  && ! useXplorNames)
	|| ((_resFlip[_resType].flags & USENEWNAMES) && (useOldNames || useXplorNames))) ) {
          validateMemo();
          _fromO = _resFlip[_resType].fromO;
          _numO  = _resFlip[_resType].numO;
          break; // found it
      }
   }
}

void FlipMemo::extendOrientations(bool val) {
   if (valid()) {
      _extendO = val;
      _fromO = _resFlip[_resType].fromO;
      _numO  = (_extendO ? _resFlip[_resType].numXO
                         : _resFlip[_resType].numO );
   }
}

bool FlipMemo::isDefaultO(int oi, SearchStrategy ss) const {
   return (ss==Mover::LOW_RES) ?
             (oi == 0) : (orientation(Mover::LOW_RES) == 0);
}

void FlipMemo::limitOrientations(bool doflip, SearchStrategy ss) {
   if (ss==Mover::LOW_RES) {

      const int halfX = ( _resFlip[_resType].numXO
			- _resFlip[_resType].numO ) / 2;
      if (! valid()) { /* don't try if not valid */ }
      else if (_extendO && (halfX > 0)) {
	 _numO = halfX;
	 _fromO = _resFlip[_resType].fromO
	        + _resFlip[_resType].numO + (doflip ? _numO : 0);
      }
      else { // most common case (HIS)
	 _numO  = _resFlip[_resType].numO / 2;
	 _fromO = _resFlip[_resType].fromO + (doflip ? _numO : 0);
      }
   }
}

void FlipMemo::finalize(int, bool, bool, bool, AtomPositions &, DotSphManager&) {
	if (isComplete()) {
		// we copy each atom and save the old locations

		for(int pi=1; pi <= _resFlip[_resType].numPnts; pi++) {
            std::map<std::string, PDBrec*>::iterator iter = _resAtoms.find(_pointName[_resType][pi]);
            
            bool found;
			if (iter != _resAtoms.end()) {
				_wrkAtom[pi] = *(iter->second);
				found = TRUE;
			}
			else {
                found = FALSE;
			}
//			bool found = _resAtoms.get(_pointName[_resType][pi], _wrkAtom[pi]);

			if (found) { _origLoc[pi] = _wrkAtom[pi].loc(); }
			else { // neighboring atom could not be found (dummy out the slot)
				_wrkAtom[pi].invalidateRecord();
				_origLoc[pi].x(-999.9).y(-999.9).z(-999.9);
			}
		}

		// finally, we generate the alternate positions for flipped H atoms
        double cur_dist;
		for(int ppi=_resFlip[_resType].fromPP;
        ppi < _resFlip[_resType].fromPP+_resFlip[_resType].numPP; ppi++) {

            if (UseNuclearDistances) {
              cur_dist = _protLoc[ppi].dist_nuclear;
            }
            else {
              cur_dist = _protLoc[ppi].dist_ecloud;
            }
            
			_origLoc[_protLoc[ppi].anum] = atomPlacementPlan::calcLoc(
				_protLoc[ppi].type,         _origLoc[_protLoc[ppi].c1],
				_origLoc[_protLoc[ppi].c2], _origLoc[_protLoc[ppi].c3],
				_origLoc[0], // *place holder* not referenced for types 3&4
				cur_dist, _protLoc[ppi].ang, _protLoc[ppi].dh);
		}
	}
}

int FlipMemo::makebumpers(std::multimap<LocBlk, BumperPoint*>& bblks, int rn, float& maxVDWrad) {
	int i = 0, an = 0;
	BumperPoint* bp;
	if (_isComplete) {
		for (i = 0; i < _resFlip[_resType].numBmpr; i++) { // regular
			const int f1 = i + 1;
			bp = new BumperPoint(_origLoc[f1], rn, an++, _wrkAtom[f1].vdwRad());
			bblks.insert(std::make_pair(LocBlk(_origLoc[f1]), bp));
//			bblks.put(LocBlk(_origLoc[f1]),
//				BumperPoint(_origLoc[f1], rn, an++, _wrkAtom[f1].vdwRad()));

			if (_wrkAtom[f1].vdwRad() > maxVDWrad) { maxVDWrad = _wrkAtom[f1].vdwRad(); }
		}
		for (i = 0; i < _resFlip[_resType].numPP; i++) { // flipped
			const int f2 = _resFlip[_resType].numPnts - i;
			bp = new BumperPoint(_origLoc[f2], rn, an++, _wrkAtom[f2].vdwRad());
			bblks.insert(std::make_pair(LocBlk(_origLoc[f2]), bp));
//			bblks.put(LocBlk(_origLoc[f2]),
//				BumperPoint(_origLoc[f2], rn, an++, _wrkAtom[f2].vdwRad()));

			if (_wrkAtom[f2].vdwRad() > maxVDWrad) { maxVDWrad = _wrkAtom[f2].vdwRad(); }
		}
	}
	return an;
}

std::list<AtomDescr> FlipMemo::getAtDescOfAllPos(float &maxVDWrad)
{
	std::list<AtomDescr> theList;
//	int i = 0;//, an = 0;
//   if (_isComplete) {
//      for (i = 0; i < _resFlip[_resType].numBmpr; i++) { // regular
//			const int f1 = i + 1;
//#ifdef DEBUGCOLLECTBUMPING
//		   cerr << "FLIPMEM: " << _origLoc[f1] << endl;
//#endif
//	 		theList->append(AtomDescr(_origLoc[f1], _wrkAtom[f1].resno(), _wrkAtom[f1].vdwRad()));
//
//	 		if (_wrkAtom[f1].vdwRad() > maxVDWrad) { maxVDWrad = _wrkAtom[f1].vdwRad(); }
//      }
//      for (i = 0; i < _resFlip[_resType].numPP; i++) { // flipped
//	 		const int f2 = _resFlip[_resType].numPnts - i;
//#ifdef DEBUGCOLLECTBUMPING
//		   cerr << "FLIPMEM: " << _origLoc[f2] << endl;
//#endif
//	 		theList->append(AtomDescr(_origLoc[f2], _wrkAtom[f2].resno(), _wrkAtom[f2].vdwRad()));
//
//	 		if (_wrkAtom[f2].vdwRad() > maxVDWrad) { maxVDWrad = _wrkAtom[f2].vdwRad(); }
//      }
//   }
   const int offO = _fromO;
   if (_resFlip[_resType].numBmpr < _resFlip[_resType].numScore)
			for (int ai = _resFlip[_resType].numScore; ai > _resFlip[_resType].numBmpr; ai--){
		     //cerr << "TEST new getAtDescOfAllPos: " << AtomDescr(_origLoc[ai], _wrkAtom[ai].resno(), _wrkAtom[ai].vdwRad()) << endl;
		     AtomDescr ad(_origLoc[ai], _wrkAtom[ai].resno(), _wrkAtom[ai].vdwRad());
		     ad.setOriginalAtomPtr( &( _wrkAtom[ ai ] ));
		     theList.push_back(ad);
		   }
   for (int j=0; j < numOrientations(Mover::LOW_RES); j++)
   {
   	for(int ai=1; ai <= _resFlip[_resType].numBmpr; ai++) {
   		if (_atomOrient[j+offO][ai] != 0) {
   			//cerr <<  "TEST new getAtDescOfAllPos: " << AtomDescr(_origLoc[_atomOrient[j+offO][ai]], _wrkAtom[ai].resno(), _wrkAtom[ai].vdwRad()) << endl;
   			AtomDescr ad(_origLoc[_atomOrient[j+offO][ai]], _wrkAtom[ai].resno(), _wrkAtom[ai].vdwRad());
   			ad.setOriginalAtomPtr( &( _wrkAtom[ ai ]));
   			theList.push_back(ad);
   		}
   	}
   }
   theList.sort();
	theList.unique();
   return theList;

}

// SJ - called from AtomPositions:setOrientation to generate the unflipped and the flipped state
bool FlipMemo::setOrientation(int oi, AtomPositions &xyz, SearchStrategy ss) {
    if (ss!=Mover::LOW_RES) { return TRUE; } // HIGH_RES uses current orientation

    // SJ - 09/10/2015 copied the static declarations to function InFlip_ModifyDAStatus
   /*static const ElementInfo * eN    = ElementInfo::StdElemTbl().element("N");
   static const ElementInfo * eNacc = ElementInfo::StdElemTbl().element("Nacc");
#ifdef AROMATICS_ACCEPT_HBONDS
   static const ElementInfo * eC    = ElementInfo::StdElemTbl().element("C");
   //static const ElementInfo * eCacc = ElementInfo::StdElemTbl().element("Car");
#endif*/

   if ((! _isComplete) || (oi < 0) || (oi >= _numO)) { return FALSE; }

    // SJ - 09/04/2015 added the flag GenerateFinalFlip so that flips are scored the same way, but finally generated using the new method
   if(!GenerateFinalFlip){
       
       Rename_Flip(oi, xyz); // SJ - 09/10/2015 moved the original code from here to this function
   }
   else if(GenerateFinalFlip){ // SJ - when GenerateFinalFlip is true - i.e. we are generating the final file and not scoring.
       bool done = FALSE;
       
       if(strcmp(_resFlip[_resType].rname,"ASN") == 0 || strcmp(_resFlip[_resType].rname,"GLN") == 0){
           if(oi == 0){ // oi ==0 means this is not the flipped orientation, got this from the array atomOrient at the top of this file
              // Rename_Flip(oi, xyz); // Nothing needs to be done, as this is not flipped, will already be in the best position when this funciton is called
               done=TRUE;
           }
           else{
              done = RotHingeDock_Flip(oi,xyz); // SJ - 09/10/2015 funciton to do the three step flip. After this funciton, the _wrkAtom will contain the new orientation.
           }
       }
       else if (strcmp(_resFlip[_resType].rname,"HIS") == 0){
           if (oi == 0 || oi == 1 || oi == 2 || oi == 6) {// SJ - 09/04/2015 these oi's correspond to unflipped state - got this from the array atomOrient at the top of this file
             //  Rename_Flip(oi, xyz); // Nothing needs to be done, as this is not flipped, will already be in the best position when this funciton is called
               done=TRUE;
           }
           else{
               done = RotHingeDock_Flip(oi,xyz); // SJ - 09/10/2015 funciton to do the three step flip. After this funciton, the _wrkAtom will contain the new orientation.
           }
       }
       
       if(!done){
           std::cerr << "Three Step Flip did not work for residue " << _resFlip[_resType].rname << " " << _wrkAtom[1].chain() << " " << _wrkAtom[1].resno() << ". Resorting to renaming flip" << std::endl;
           Rename_Flip(oi, xyz);
       }
   }
   rememberOrientation(oi);
   return TRUE;
}

// SJ - TODO explain more clearly
/*SJ - 09/10/2015 function to do the standard flip, given the correct orientation. 
 Copied the original code from setOrinetations here. 
 This function will either flip or not flip depending on the value of the orientation passed. 
 No need for explicit check for the orientation, as it is hardcoded in the atomOrient array on top of this file*/
void FlipMemo::Rename_Flip(int orientation, AtomPositions & xyz) // SJ 10/14/2015 - making this function exactly like the old setOrientation()
{
    static const ElementInfo * eN    = ElementInfo::StdElemTbl().element("N");
    static const ElementInfo * eNacc = ElementInfo::StdElemTbl().element("Nacc");
#ifdef AROMATICS_ACCEPT_HBONDS
    static const ElementInfo * eC    = ElementInfo::StdElemTbl().element("C");
    //static const ElementInfo * eCacc = ElementInfo::StdElemTbl().element("Car");
#endif

    const int offO = _fromO;
    
     /**FLIP**/
    // SJ - this is where the coordinates of the atoms are exchanged. for Hs the coordinates have already been recalcualted and stored for the flipped position in FlipMemo::Finalize();
    for(int ai=1; ai <= _resFlip[_resType].numBmpr; ai++) {
        if (_atomOrient[orientation+offO][ai] != 0) { // SJ - if this atom exists in this orientation, basically for diff protonated states of HIS
            Point3d lastloc = _wrkAtom[ai].loc();
            _wrkAtom[ai].revalidateRecord();
            _wrkAtom[ai].loc(_origLoc[_atomOrient[orientation+offO][ai]]); // SJ - atoms coordinates exchanged according to atomOrient array at the top of the file. atomOrient contains what atoms are to be swapped with what, according to the value of the oi that designates to flip or not to flip.
            xyz.reposition(lastloc,_wrkAtom[ai]);
            
            if (_wrkAtom[ai].elem().atno() == 7) {
                if (_atomDAflags[orientation+offO][ai] < 0) {
                    _wrkAtom[ai].elem(*eNacc);
                }
                else { _wrkAtom[ai].elem(*eN); }
            }
#ifdef AROMATICS_ACCEPT_HBONDS
            else if (_wrkAtom[ai].elem().atno() == 6) {
                _wrkAtom[ai].elem(*eC); // apl - 10/19/2006 - His carbons no longer aromatic
                //if (_atomDAflags[oi+offO][ai] < 0) {
                //  _wrkAtom[ai].elem(*eCacc);
                //}
                //else { _wrkAtom[ai].elem(*eC); }
            }
#endif

        }
        else { // switch off but do not rub out completely
            _wrkAtom[ai].partiallyInvalidateRecord();
        }
    }
    
    // This was done after the original code as well, so just following the same thing
   /* for(int ai=1; ai <= _resFlip[_resType].numBmpr; ai++) {
        if (_atomOrient[orientation+offO][ai] != 0) {
            // *** notice: the D/A status of nitrogen and carbon are modified here ***
            InFlip_ModifyDAStatus(ai,orientation,offO); // SJ - 09/10/2015 moved the original code from setOrientation to this function
        }
    }*/
    
    return;
}

// SJ - 09/10/2015 function to do the three step flip
bool FlipMemo::RotHingeDock_Flip(int orientation, AtomPositions & xyz)
{
    Point3d newLoc, lastloc[FMmaxAtomSlots+1];
    const int offO = _fromO;

    Point3d p1,p2,p3,p4,a1,b1,c1,a2,b2,c2,normal1,normal2,rotaxis;
    double rotangle=0;
    
    /**SETUP**/
    
    //storing the last location of all the atoms, and copying the original coordinates back into the working copy, because you have to flip from the original coordinates
    for(int pi=1; pi <= _resFlip[_resType].numPnts; pi++) {
        lastloc[pi] = _wrkAtom[pi].loc();
        _wrkAtom[pi].loc(_origLoc[pi]);
    }
    
    // revalidate the records for the atoms that are in this orientation, and invalidate the ones that are not.
    // basically for diff protonated states of HIS - taken directly from the standard flip code
    for(int ai=1; ai <= _resFlip[_resType].numBmpr; ai++) {
        if (_atomOrient[orientation+offO][ai] != 0) {
            _wrkAtom[ai].revalidateRecord();
        }
        else { // switch off but do not rub out completely
            _wrkAtom[ai].partiallyInvalidateRecord();
        }
    }
    
    /**ROTATION**/
    
    //getting the axis for the first 180 degree rotation
    // indices from _dockAtomIndex array at the top of this file
    p1 = _origLoc[_dockAtomIndex[_resType][1]]; //CG for GLN, CB for ASN, HIS
    p2 = _origLoc[_dockAtomIndex[_resType][2]]; //CD for GLN, CG for ASN, HIS
    
    //rotate the atoms involved in the 180 degree roation
    for(int ai=1; ai <= _resFlip[_resType].numBmpr; ai++) {
        if (_atomOrient[orientation+offO][ai] != 0) {
            
            newLoc = _origLoc[ai].rotate(180,p1,p2); // rotate the original position of the atom
            /*if(newLoc == NULL) // rotation did not work, return FALSE, TODO: need to modify the criteria for failure
                return FALSE;*/
            _wrkAtom[ai].loc(newLoc);
       }
    }
    
    /**HINGE**/
    
    //getting the points that define the plane for the terminal groups
    //original coordinates
    a1 = _origLoc[_dockAtomIndex[_resType][2]]; // CD for GLN, CG for ASN, HIS
    b1 = _origLoc[_dockAtomIndex[_resType][3]]; // OE1 for GLN, OD1 for ASN, CE1 for HIS
    c1 = _origLoc[_dockAtomIndex[_resType][4]]; // NE2 for GLN, HIS, ND2 for ASN
    //new coordinates
    a2 = _wrkAtom[_dockAtomIndex[_resType][2]].loc();
    b2 = _wrkAtom[_dockAtomIndex[_resType][3]].loc();
    c2 = _wrkAtom[_dockAtomIndex[_resType][4]].loc();
    
    //computing normal to the plane by cross product of the two vectors in the plane
    normal1=cross(makeVec(b1, a1), makeVec(c1, a1)); // normal to original plane
    normal2=cross(makeVec(b2, a2), makeVec(c2, a2)); // normal to new plane
    normal1.normalize();
    normal2.normalize();
    
    //rotation axis is the cross product of two normals, angle is acos of the dot product
    rotaxis=cross(normal1,normal2);
    rotangle=acos(dot(normal1,normal2))*180/PI; // converting to degrees
    rotangle=180-rotangle; // as the normals are facing each other, this needs to be corrected.
    
    //rotate the new coordinates of the atoms involved in the hinge about the rotaxis
    for(int ai=1; ai <= _resFlip[_resType].numBmpr; ai++) {
        if (_atomOrient[orientation+offO][ai] != 0) { 
            newLoc = _wrkAtom[ai].loc().rotate(rotangle,a1,a1+rotaxis); // as a1 is always on the line of intersection of the plane, i.e. rotation axis vector
            //if(newLoc == NULL) // rotation did not work, return FALSE, TODO: need to modify the criteria for failure
            // return FALSE;
            _wrkAtom[ai].loc(newLoc);
        }
    }
    
    /**DOCK**/
    
    // first rotation - getting the points for the two vectors
    p2 = _origLoc[_dockAtomIndex[_resType][5]]; // CA
    p3 = _origLoc[_dockAtomIndex[_resType][4]]; // original ND2 for ASN, NE2 for GLN and HIS
    p4 = _wrkAtom[_dockAtomIndex[_resType][3]].loc(); // mobile OD1 for ASN, OE1 for GLN, CE1 for HIS
    
    //axis is cross product of the vectors, angle is acos of dot product
    rotaxis=cross(makeVec(p2,p3),makeVec(p2,p4));
    rotangle=acos(dot(makeVec(p2,p3),makeVec(p2,p4)))*180/PI; // makeVec returns normalized vectors
    rotangle=360-rotangle; // rotation has to be clockwise around the axis
    
    //rotating all the atoms that are included in the CA dock
    //_dockAtomIndex[_resType][0] contains the number of atoms to rotate
    
    for(int i=1; i <= _dockAtomIndex[_resType][0]; i++){
        
        if(_wrkAtom[_dockAtomIndex[_resType][i]].valid()){ // if this atom has not be partiallyInvalidated for this orientation
            newLoc = _wrkAtom[_dockAtomIndex[_resType][i]].loc().rotate(rotangle,p2,p2+rotaxis); // since the rotation has to be centered around the rotaxis passing through the CA
            //if(newLoc == NULL) // rotation did not work, return FALSE, TODO: need to modify the criteria for failure
            // return FALSE;
            _wrkAtom[_dockAtomIndex[_resType][i]].loc(newLoc);
        }
    }
    
    //second rotation - getting the four points for the dihedral angle, p2 and p3 remains the same as above
    p1=_wrkAtom[_dockAtomIndex[_resType][4]].loc(); // mobile ND2 for ASN, NE2 for GLN, HIS
    p4=_origLoc[_dockAtomIndex[_resType][3]]; // original OD1 for ASN, OE1 for GLN, CE1 for HIS
    
    //angle is equal to the dihedral angle, rotation axis is p2->p3
    rotangle=dihedral(p1,p2,p3,p4);
    
    //rotating all the atoms that are included in the CA dock
    //_dockAtomIndex[_resType][0] contains the number of atoms to rotate
    for(int i=1; i <= _dockAtomIndex[_resType][0]; i++){
        
        if(_wrkAtom[_dockAtomIndex[_resType][i]].valid()){ // if this atom has not be partiallyInvalidated for this orientation
        
            newLoc = _wrkAtom[_dockAtomIndex[_resType][i]].loc().rotate(rotangle,p2,p3);
            //if(newLoc == NULL) // rotation did not work, return FALSE, TODO: need to modify the criteria for failure
            // return FALSE;
            _wrkAtom[_dockAtomIndex[_resType][i]].loc(newLoc);
        }
    }
    
    /**FINAL STEPS**/ // These were done after the standard flip in the original code as well, so just following the same thing
    
    //reposition all the atoms for which the coordinates have changed
    for(int i=1; i <= _dockAtomIndex[_resType][0]; i++){
        
        if(_wrkAtom[_dockAtomIndex[_resType][i]].valid()){ // if this atom has not be partiallyInvalidated for this orientation
            xyz.reposition(lastloc[_dockAtomIndex[_resType][i]], _wrkAtom[_dockAtomIndex[_resType][i]]);
        }
    }
    
    for(int ai=1; ai <= _resFlip[_resType].numBmpr; ai++) {
        if (_atomOrient[orientation+offO][ai] != 0) {
            // *** notice: the D/A status of nitrogen and carbon are modified here ***
            InFlip_ModifyDAStatus(ai,orientation,offO); // SJ - 09/10/2015 moved the original code from setOrientation to this function
        }
    }
    
    return TRUE;
}

// SJ - 09/10/2015 function to modify the donor acceptor status of N and C atoms after flips have happened. Copied the original code from setOrientations and made a function because it has to called from multiple places.
void FlipMemo::InFlip_ModifyDAStatus(int ai, int orientation, int offO){
    
    static const ElementInfo * eN    = ElementInfo::StdElemTbl().element("N");
    static const ElementInfo * eNacc = ElementInfo::StdElemTbl().element("Nacc");
#ifdef AROMATICS_ACCEPT_HBONDS
    static const ElementInfo * eC    = ElementInfo::StdElemTbl().element("C");
     //static const ElementInfo * eCacc = ElementInfo::StdElemTbl().element("Car");
#endif
    
    if (_wrkAtom[ai].elem().atno() == 7) {
        if (_atomDAflags[orientation+offO][ai] < 0) {
            _wrkAtom[ai].elem(*eNacc);
        }
        else { _wrkAtom[ai].elem(*eN); }
    }
#ifdef AROMATICS_ACCEPT_HBONDS
    else if (_wrkAtom[ai].elem().atno() == 6) {
        _wrkAtom[ai].elem(*eC); // apl - 10/19/2006 - His carbons no longer aromatic
        //if (_atomDAflags[oi+offO][ai] < 0) {
        //  _wrkAtom[ai].elem(*eCacc);
        //}
        //else { _wrkAtom[ai].elem(*eC); }
    }
#endif
    return;
}

std::string FlipMemo::describeOrientation() const {
   const int oi = orientation();
   const char *dscr =
       ((! valid()) || (oi < 0) || (oi >= _numO))
      ? "" : _orientDescr[oi + _fromO];
   return dscr;
}

int FlipMemo::flipState() const {
   const int oi = orientation();
   return ((! valid()) || (oi < 0) || (oi >= _numO))
      ? 0 : _oFlipState[oi + _fromO];
}

bool FlipMemo::markFlipAtoms() {
   if (_wrkAtom[1].atomNameModified() || (flipState() == 0)) {
      return FALSE;
   }
   for(int ai=1; ai <= _resFlip[_resType].numHeavy; ai++) {
      _wrkAtom[ai].annotation("flip");
   }
   return TRUE;
}

double FlipMemo::determineScore(AtomPositions &xyz,	DotSphManager& dotBucket,
								int nBondCutoff, float probeRadius, float pmag,
								double& penalty, float &bumpScore, float &hbScore,
								bool& hasBadBump) {

	bumpScore  = 0.0;
	hbScore    = 0.0;
	hasBadBump = FALSE;

	const double maxVDWrad = ElementInfo::StdElemTbl().maxExplicitRadius();
	const int numA = numScoreAtoms();

	double scoreThisO = 0.0;
	for(int j=0; j < numA; j++) {
		PDBrec thisAtom;
		float bumpSubScore = 0.0;
		float hbSubScore   = 0.0;
		bool  subBadBump   = FALSE;

		if (locAtom(j, thisAtom)) {
			std::list<PDBrec*> bnded = neighbors(j, nBondCutoff);

			//std::cerr << "Scoring atom: " << j << " " << thisAtom.getAtomDescr() << std::endl;
			//for (std::list< PDBrec* >::iterator iter = bnded.begin(); iter != bnded.end(); ++iter )
			//{
			//	std::cerr << "bonded for FlipMemo: " << nBondCutoff << " " << (*iter)->getAtomDescr() << std::endl;
			//}
            
			double val = xyz.atomScore(thisAtom, thisAtom.loc(),
				thisAtom.vdwRad() + probeRadius + maxVDWrad,
				//apl procrastinate nearby list computation until AtomPositions decides to score
				bnded,
				dotBucket.fetch(thisAtom.vdwRad()),
				probeRadius, FALSE,
				bumpSubScore, hbSubScore, subBadBump);
            
///////////////////////////////////////////////////////////////////////////////////////////
// Delete bnded!
//			std::for_each(bnded.begin(), bnded.end(), DeleteObject());
///////////////////////////////////////////////////////////////////////////////////////////

			bumpScore  += bumpSubScore;
			hbScore    += hbSubScore;
			if (subBadBump) { hasBadBump = TRUE; }

#ifdef DEBUGSUBSCORE
            std::cerr << "\t:" << thisAtom.recName()
				<<":"<< describeOrientation()
				<< ": bump=" << bumpSubScore
				<< ", HB=" << hbSubScore
				<< ((subBadBump==TRUE)?", BADBUMP":"")
				//      << ", tot=" << val
				<< "\t"
				//      << thisAtom.loc()
            << std::endl;
            
#endif
			scoreThisO += val;
		}
	}
	penalty = orientationPenalty(pmag);
	initializeScoreIfNotSet(scoreThisO, hasBadBump);

	return scoreThisO;
}

void FlipMemo::insertAtom(PDBrec* r) {
   if (valid()) {
	   _resAtoms.insert(std::make_pair(r->atomname(), r));
      fillAtomAndLocVectors(); // try and find required atoms
   }
}

// externally, we can refer to the atom num in the range [0..numScore)
//SJ - this function puts the PDBRec for the specific atom (specified by na) in the outrec
bool FlipMemo::locAtom(int na, PDBrec& outrec) const {
   bool rc = FALSE;
   if (_isComplete) {
      if (na >= 0 || na < _resFlip[_resType].numScore) {
	 outrec = _wrkAtom[na+1];
	 if (outrec.valid() && (! outrec.hasProp(IGNORE)) ) {
	    rc = TRUE;  // mark the invalid
	 }
      }
   }
   return rc;
}

// externally, na is the atom num in the range [0..numScore)
std::list<PDBrec*> FlipMemo::neighbors(int na, int nbdist) const {
	std::list<PDBrec*> nbhd;
	int nnr = -1;
	if (_isComplete) { // don't bother if incomplete
		for(int i=_resFlip[_resType].from3B;
		i < _resFlip[_resType].from3B+_resFlip[_resType].numScore; i++) {

			if (_bondedLimits[i].anum == na+1) { nnr = i; break; }
		}
	}
	if (nnr >= 0) { // found the neighbor record
		int nbcnt = 0;
		switch(nbdist) {
		case 1:  nbcnt = _bondedLimits[nnr].b1; break;
		case 2:  nbcnt = _bondedLimits[nnr].b2; break;
		default: nbcnt = _bondedLimits[nnr].b3; break;
		}
		for(int j=0; j < nbcnt; j++) {
			const int k = _bondedSet[nnr][j];
			if (k <= 0
				|| k > _resFlip[_resType].numPnts) { break; } // error...

			if (_wrkAtom[k].valid() && (! _wrkAtom[k].hasProp(IGNORE)) ) {
				PDBrec* temp = new PDBrec(_wrkAtom[k]);
				nbhd.push_front(temp);
			}
		}
	}
	return nbhd;
}

// SJ - this function checks if this residue is the type to be checked for flips, and puts the altcodes for the atoms to be involved in the argument sch. The type of residues to be flipped, and the involved atoms are stored in the class FlipMemo. Look at the top of this file.
void FlipMemo::altCodes(const ResBlk& rblk,	bool useXplorNames, bool useOldNames, bool bbModel,
						std::list<char>& sch) {
	char ch, buf[10];
	int rt = 0, k = 0, cursor = 0;
	bool isInResidueSet = FALSE, dupalt = FALSE;

	const char *resname = rblk.firstRec().resname();

	for(rt = 0; _resFlip[rt].rname; rt++) { // SJ checks if the residue name is ASN, GLN, or HIS
        
	        if ( (strcmp(_resFlip[rt].rname, resname) == 0)
	      		&& !(((_resFlip[rt].flags & USEOLDNAMES) && ! useOldNames)
	        	|| ((_resFlip[rt].flags & XPLORNAME)  && ! useXplorNames)
	        	|| ((_resFlip[rt].flags & USENEWNAMES) && (useOldNames || useXplorNames))) ) {
			isInResidueSet = TRUE;
			break; // residue type is one of those we are concerned with
		}
	}
	if (isInResidueSet) { // SJ - if the residue is ASN, GLN, or HIS
		std::multimap<std::string, PDBrec*> pdb = rblk.atomIt();
		std::string key;
		std::multimap<std::string, PDBrec*>::const_iterator pdbit = pdb.begin();
		PDBrec* atsq = NULL;
		while(pdbit != pdb.end()) {
			key = pdbit->first;
			for (; pdbit != pdb.end() && pdbit->first == key; ++pdbit) {
				atsq = pdbit->second;
				bool foundname = FALSE;
				for(int i=_resFlip[rt].fromScat;
				i < _resFlip[rt].fromScat+_resFlip[rt].numScat; i++) { // SJ - for each atom that takes part in the flip
                    if (strcmp(_pointName[rt][i], atsq->atomname()) == 0) {
                        foundname = TRUE;
						break;
					}
				}
				if (foundname) { // SJ - look for alternate conformations code
					ch = toupper(atsq->alt());
					if (ch != ' ') {
						dupalt = FALSE;
						for(k = 0; k < cursor; k++) {
							if (ch == buf[k]) { dupalt = TRUE; break; }
						}
						if (! dupalt) {
							buf[cursor++] = ch;
						}
					}
				}
			}
		}
		if (cursor < 1) { sch.push_front(' '); } // no alt codes
		else {
			for(k = 0; k < cursor; k++) { // at least one alt code
				sch.push_front(buf[k]);
			}
		}
	}
}

//Is this a heavy atom which is an HBond donor or acceptor when flipped?

bool FlipMemo::isHBDonorOrAcceptorFlipped(const PDBrec& a, bool useXplorNames, bool useOldNames, bool bbModel) {
   int rt = 0;
   bool isFlipableResidue = FALSE;

   const char *resname = a.resname();

   for(rt = 0; _resFlip[rt].rname; rt++) {
      if ( (strcmp(_resFlip[rt].rname, resname) == 0)
              && !(((_resFlip[rt].flags & USEOLDNAMES) && ! useOldNames)
              || ((_resFlip[rt].flags & XPLORNAME)  && ! useXplorNames)
              || ((_resFlip[rt].flags & USENEWNAMES) && (useOldNames || useXplorNames))) ) {
	 isFlipableResidue = TRUE; break;
      }
   }
   if (isFlipableResidue) {
      for(int i=_resFlip[rt].fromScat;
	 i < _resFlip[rt].fromScat+_resFlip[rt].numHeavy; i++) {
	 if (strcmp(_pointName[rt][i], a.atomname()) == 0) {
	    return TRUE;
	 }
      }
   }
   return FALSE;
}

void FlipMemo::dropBondedFromBumpingListForPDBrec(
		std::list< PDBrec * > & bumping,
		PDBrec* atom,
		int nBondCutoff
) const
{
	//std::cerr << "Dropping bonded from bumping list for " << atom->getAtomDescr() << std::endl;
	//for ( std::list< PDBrec * >::iterator iter = bumping.begin();
	//	iter != bumping.end(); ++iter )
	//{
	//	std::cerr << "Bumping: " << nBondCutoff << " " << (*iter)->getAtomDescr() << std::endl;
	//}

	int atom_number = findAtom( atom );
	//std::cerr << "Atom number: " << atom_number << std::endl;

	std::list<PDBrec*> bnded = neighbors(atom_number, nBondCutoff);
	for (std::list< PDBrec* >::iterator iter = bumping.begin();
		iter != bumping.end(); )
	{
		std::list< PDBrec* >::iterator iter_next = iter;
		++iter_next;
		//if ( std::find( bnded.begin(), bnded.end(), *iter ) != bnded.end() )
		//{
		//	bumping.erase( iter );
		//}
		for (std::list< PDBrec* >::const_iterator constiter = bnded.begin();
			constiter != bnded.end(); ++constiter)
		{
			//std::cerr << "Compairing against: " << (*constiter)->getAtomDescr() << std::endl;

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

int FlipMemo::findAtom( PDBrec * atom ) const
{
	for (int ii = 1; ii <= FMmaxAtomSlots; ++ii)
	{
		if ( &_wrkAtom[ ii ] == atom )
			return ii - 1; //FIXING OFF BY ONE ERROR
	}
	std::cerr << "Critical error in FlipMemo::findAtom( " << atom << " ).  Atom not found." << std::endl;
	for (int ii = 1; ii <= FMmaxAtomSlots; ++ii)
	{
		std::cerr << " &_wrkAtom[ ii ] " << &_wrkAtom[ ii ]  << std::endl;
	}
	exit(2);
        return 0; // to avoid warnings
}

double FlipMemo::orientationPenalty(float pmag) const {
   if (_wrkAtom[1].atomNameModified()) { pmag = 0.0; }
   const int oi = orientation();
   const double pf = ((! valid()) || (oi < 0) || (oi >= _numO))
      ? 0.0
      : _orientationPenalty[oi + _fromO];
   return -pmag * pf;
}

void FlipMemo::fillAtomAndLocVectors() {
	if (!valid()) { return; }

	_wrkAtom[0].invalidateRecord(); // dummy slot zero
	_origLoc[0].x(-999.9).y(-999.9).z(-999.9);

	// first a scan to see if the required atoms have been gathered
	// (i.e., those required to calc the position flipped protons)

	for(int i=_resFlip[_resType].fromScat;
	i < _resFlip[_resType].fromScat+_resFlip[_resType].numScat; i++) {
		if (_resAtoms.find(_pointName[_resType][i]) == _resAtoms.end()) {
			_isComplete = FALSE;
			return;
		}

	}
	_isComplete = TRUE; // good enough
}
