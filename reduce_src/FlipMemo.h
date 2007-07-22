// name: FlipMemo.h
// author: J. Michael Word
// date written: 2/7/98
// purpose: Interface for FlipMemo

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

#ifndef FLIPMEMO_H
#define FLIPMEMO_H 1

#ifdef OLD_STD_HDRS
#include <limits.h>
#else
#include <climits>
#endif

#include "PDBrec.h"
#include "ResBlk.h"
#include "DotSph.h"
#include "BumperPoint.h"
#include "Mover.h"
#include "neighbors.h"
#include "utility.h"

struct ResFlipInfo_t {
   const char* rname;  // residue name
   int fromO;    // first orientation
   int numO;     // number of orientations
   int numXO;    // number of orient. w/extended
   int fromPP;   // first proton placement record
   int numPP;    // number of proton placement recs
   int from3B;   // first 1-3 bond neighbor record
   int numScore; // number of atoms to be scored
   int numBmpr;  // num orig atoms which bump or re-orient
   int numHeavy; // number of heavy atoms moved
   int fromScat; // first scorable atom
   int numScat;  // number of scorable atoms
   int numPnts;  // number
   int flags;
};
struct ProtonLoc_t {
   const char* rname;
   const char* aname;
   int type, anum, c1, c2, c3;
   float dist, ang, dh;
};
struct BondLimits_t {
   int anum, b1, b2, b3;
};

// flip Memo constants used to size data tables
const int FMmaxAtomSlots     =20;
const int FMmaxBondedAtoms   =16;
const int FnumScoreableAtoms =19;
const int FMnumFlipOrient    =10;
const int FMnumNewLocs       = 8;
const int FMnumFlipRes       = 8;

class FlipMemo: public Mover {
public:
   FlipMemo(const char *resname, bool useXplorNames, bool useOldNames, bool bbModel);
   virtual ~FlipMemo() {
   }

   virtual Mover::MemoType type() { return Mover::FLIP; }
   virtual bool isComplete() const { return _isComplete; }
   virtual bool hasHires() const { return FALSE; }
   virtual bool canFlip() const { return TRUE; }
   virtual bool canRotate() const { return FALSE; }
   virtual int flipState() const;
   virtual bool markFlipAtoms();
   virtual void finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel,
                         AtomPositions &xyz, DotSphManager& dotBucket);
   virtual int makebumpers(std::multimap<LocBlk, BumperPoint*>& bbins,
                           int n, float& maxVDWrad);
   virtual std::list<AtomDescr> getAtDescOfAllPos(float &maxVDWrad);                        
   virtual const PDBrec& exampleAtom() const { return _wrkAtom[1]; }

   virtual int numOrientations(
	 Mover::SearchStrategy ss=Mover::LOW_RES) const {
      return valid() ? ((ss==Mover::LOW_RES) ? _numO : 1) : 0;
   }
   virtual void limitOrientations(bool doflip,
	 SearchStrategy ss=Mover::LOW_RES);
   virtual bool setOrientation(int oi, AtomPositions &xyz,
	 SearchStrategy ss=Mover::LOW_RES);
   virtual bool isDefaultO(int oi,
	 SearchStrategy ss=Mover::LOW_RES) const;
   virtual std::string describeOrientation() const;
   virtual double determineScore(AtomPositions &xyz,
      DotSphManager& dotBucket, int nBondCutoff,
      float probeRadius, float pmag, double& penalty,
      float &bumpScore, float &hbScore, bool& hasBadBump);

   void extendOrientations(bool val);

   void insertAtom(PDBrec* r);

   int numScoreAtoms() const {
      return valid() ? _resFlip[_resType].numScore : 0;
   }
   bool locAtom(int ai, PDBrec& outrec) const;

   std::list<PDBrec*> neighbors(int na, int nbdist) const;

   static void altCodes(const ResBlk& rblk, bool useXplorNames, bool useOldNames, bool bbModel, std::list<char>& sch);
   static bool isHBDonorOrAcceptorFlipped(const PDBrec& a, bool useXplorNames, bool useOldNames, bool bbModel);

   virtual void setHydAngle(double, AtomPositions &) {/*do nothing*/}
   
   virtual void dropBondedFromBumpingListForPDBrec( std::list< PDBrec * > & bumping, PDBrec* atom, int nBondCutoff  ) const;
private:
   double orientationPenalty(float pmag) const;
   void fillAtomAndLocVectors(); // utility routines
   int findAtom( PDBrec* atom ) const; // which atom in the list is a particular PDBrec?

   FlipMemo(const FlipMemo& m); // assign and copy are not implemented
   FlipMemo& operator=(const FlipMemo& m);

   std::map<std::string, PDBrec*> _resAtoms;
   int         _resType;
   bool        _isComplete;   // have all the atoms been gathered?
   bool        _extendO;      // extend the number of orientations?

   int         _fromO;        // first orientation
   int         _numO;         // number of orientations?

   PDBrec  _wrkAtom[FMmaxAtomSlots+1]; // the atoms we will work with
   Point3d _origLoc[FMmaxAtomSlots+1]; // original xyz locations of atoms

   // tables which describe how to flip particular residues
   static const char*         _pointName[FMnumFlipRes+1][FMmaxAtomSlots+1];
   static const ResFlipInfo_t _resFlip[FMnumFlipRes+1];
   static const ProtonLoc_t   _protLoc[FMnumNewLocs+1];
   static const BondLimits_t  _bondedLimits[FnumScoreableAtoms];
   static const int           _bondedSet[FnumScoreableAtoms][FMmaxBondedAtoms];
   static const int           _atomOrient[FMnumFlipOrient][FMmaxAtomSlots+1];
   static const int           _atomDAflags[FMnumFlipOrient][FMmaxAtomSlots+1];
   static const char*         _orientDescr[FMnumFlipOrient];
   static const double        _orientationPenalty[FMnumFlipOrient];
   static const int           _oFlipState[FMnumFlipOrient];
};
#endif
