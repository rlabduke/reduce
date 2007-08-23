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

const ResFlipInfo_t FlipMemo::_resFlip[FMnumFlipRes+1] = {
// --- rname, fromO,numO,numXO, fromPP,numPP, from3B,numScore,numBmpr,numHeavy,
// ---        fromScat,numScat, numPnts, flags

#ifdef AROMATICS_ACCEPT_HBONDS
 {"HIS", 0,6,8, 0,4, 0,9,9,4, 1,9, 20,           0},
 {"HIS", 0,6,8, 0,4, 0,9,9,4, 1,9, 20, USENEWNAMES},
#else
 {"HIS", 0,6,8, 0,4, 0,8,8,4, 1,9, 20,           0},
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
 {"HIS", " HD1", 4, 17, 2, 9, 4, 1.0,   0.0,   0.0},
 {"HIS", " HD2", 4, 18, 1, 3, 9, 1.1,   0.0,   0.0},
 {"HIS", " HE1", 4, 19, 4, 2, 3, 1.1,   0.0,   0.0},
 {"HIS", " HE2", 4, 20, 3, 4, 1, 1.0,   0.0,   0.0},
 {"ASN", "1HD2", 3, 13, 1, 5, 2, 1.0, 120.0,   0.0},
 {"ASN", "2HD2", 3, 14, 1, 5, 2, 1.0, 120.0, 180.0},
 {"GLN", "1HE2", 3, 13, 1, 5, 2, 1.0, 120.0,   0.0},
 {"GLN", "2HE2", 3, 14, 1, 5, 2, 1.0, 120.0, 180.0},
 { NULL,   NULL, 0,  0, 0, 0, 0, 0.0,   0.0,   0.0}
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
                               [FMmaxAtomSlots+1] = {
#ifdef AROMATICS_ACCEPT_HBONDS
 {0, 1, 2, 3, 4, 0, 6, 7, 8, 9}, // HIS...
 {1, 1, 2, 3, 4, 5, 6, 7, 0, 9},
 {2, 1, 2, 3, 4, 5, 6, 7, 8, 9},
 {3, 2, 1, 4, 3, 0,18,19,20, 9},
 {4, 2, 1, 4, 3,17,18,19, 0, 9},
 {5, 2, 1, 4, 3,17,18,19,20, 9},
 {6, 1, 2, 3, 4, 0, 6, 7, 0, 9}, // rarely used double deprotonated states
 {7, 2, 1, 4, 3, 0,18,19, 0, 9},
#else
 {0, 1, 2, 3, 4, 0, 6, 7, 8}, // HIS...
 {1, 1, 2, 3, 4, 5, 6, 7, 0},
 {2, 1, 2, 3, 4, 5, 6, 7, 8},
 {3, 2, 1, 4, 3, 0,18,19,20},
 {4, 2, 1, 4, 3,17,18,19, 0},
 {5, 2, 1, 4, 3,17,18,19,20},
 {6, 1, 2, 3, 4, 0, 6, 7, 0}, // rarely used double deprotonated states
 {7, 2, 1, 4, 3, 0,18,19, 0},
#endif
 {0, 1, 2, 3, 4, 5},          // ASN, GLN...
 {1, 2, 1,13,14, 5}
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

		for(int ppi=_resFlip[_resType].fromPP;
        ppi < _resFlip[_resType].fromPP+_resFlip[_resType].numPP; ppi++) {
			
			_origLoc[_protLoc[ppi].anum] = atomPlacementPlan::calcLoc(
				_protLoc[ppi].type,         _origLoc[_protLoc[ppi].c1],
				_origLoc[_protLoc[ppi].c2], _origLoc[_protLoc[ppi].c3],
				_origLoc[0], // *place holder* not referenced for types 3&4
				_protLoc[ppi].dist, _protLoc[ppi].ang, _protLoc[ppi].dh);
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

bool FlipMemo::setOrientation(int oi, AtomPositions &xyz, SearchStrategy ss) {
   if (ss!=Mover::LOW_RES) { return TRUE; } // HIGH_RES uses current orientation

   static const ElementInfo * eN    = ElementInfo::StdElemTbl().element("N");
   static const ElementInfo * eNacc = ElementInfo::StdElemTbl().element("Nacc");
#ifdef AROMATICS_ACCEPT_HBONDS
   static const ElementInfo * eC    = ElementInfo::StdElemTbl().element("C");
   //static const ElementInfo * eCacc = ElementInfo::StdElemTbl().element("Car");
#endif

   if ((! _isComplete) || (oi < 0) || (oi >= _numO)) { return FALSE; }

   const int offO = _fromO;

   for(int ai=1; ai <= _resFlip[_resType].numBmpr; ai++) {
      if (_atomOrient[oi+offO][ai] != 0) {
	 Point3d lastloc = _wrkAtom[ai].loc();
	 _wrkAtom[ai].revalidateRecord();
	 _wrkAtom[ai].loc(_origLoc[_atomOrient[oi+offO][ai]]);
	 xyz.reposition(lastloc, _wrkAtom[ai]);

// *** notice: the D/A status of nitrogen and carbon are modified here ***
	 if (_wrkAtom[ai].elem().atno() == 7) {
	    if (_atomDAflags[oi+offO][ai] < 0) {
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
   rememberOrientation(oi);
   return TRUE;
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
			cerr << "\t:" << thisAtom.recName()
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

void FlipMemo::altCodes(const ResBlk& rblk,	bool useXplorNames, bool useOldNames, bool bbModel, 
						std::list<char>& sch) {
	char ch, buf[10];
	int rt = 0, k = 0, cursor = 0;
	bool isInResidueSet = FALSE, dupalt = FALSE;

	const char *resname = rblk.firstRec().resname();

	for(rt = 0; _resFlip[rt].rname; rt++) {
	        if ( (strcmp(_resFlip[rt].rname, resname) == 0)
	      		&& !(((_resFlip[rt].flags & USEOLDNAMES) && ! useOldNames)
	        	|| ((_resFlip[rt].flags & XPLORNAME)  && ! useXplorNames)
	        	|| ((_resFlip[rt].flags & USENEWNAMES) && (useOldNames || useXplorNames))) ) {
			isInResidueSet = TRUE;
			break; // residue type is one of those we are concerned with
		}
	}
	if (isInResidueSet) {
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
				i < _resFlip[rt].fromScat+_resFlip[rt].numScat; i++) {

					if (strcmp(_pointName[rt][i], atsq->atomname()) == 0) {
						foundname = TRUE;
						break;
					}
				}
				if (foundname) {
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
	exit(1);
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
