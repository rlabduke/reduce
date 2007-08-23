// name: Mover.h
// author: J. Michael Word
// date written: 2/7/97
// purpose: Interface for Mover

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
#pragma warning(disable:4800) 
#endif

#ifndef MOVER_H
#define MOVER_H 1

#ifdef OLD_STD_HDRS
#include <limits.h>
#else
#include <climits>
#endif

#include <map>
#include <functional>
#include "PDBrec.h"
#include "DotSph.h"
#include "BumperPoint.h"
#include "neighbors.h"
#include "utility.h"


class AtomPositions;
class Mover;
typedef Mover* MoverPtr;

const double LowestMoverScore = -9.9E99;

// (abstract) base class for movement memos
class Mover {
protected:
   Mover() :  _bestScore(LowestMoverScore),
	      _initScore(LowestMoverScore),
	      _bestHasBadBump(FALSE), _initHasBadBump(FALSE),
              _initIsSet(FALSE), _ok(FALSE), _adj(TRUE) {
      _orientation[0] = 0;
      _orientation[1] = 0;
      _bestOrientation[0] = 0;
      _bestOrientation[1] = 0;
      resetFlipMaxScore(0);
      resetFlipMaxScore(1);
   }
public:
   virtual ~Mover() {
//	   for (std::map<std::string, MoverPtr>::const_iterator i = _links.begin(); i != _links.end(); ++i)
//		   delete i->second;
   }

   enum MemoType { ROTATE_METHYL, ROTATE_DONOR, FLIP };
   enum SearchStrategy { LOW_RES, HIGH_RES };

   virtual MemoType type() = 0;
   virtual bool isComplete() const = 0;
   virtual bool hasHires() const = 0;
   virtual bool canFlip() const = 0;
   virtual bool canRotate() const = 0;
   virtual int flipState() const = 0;
   virtual bool markFlipAtoms() = 0;
   virtual void finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel, AtomPositions &xyz, DotSphManager& dotBucket) = 0;
   virtual const PDBrec& exampleAtom() const = 0;
   virtual int makebumpers(std::multimap<LocBlk, BumperPoint*>& bbins,
                           int n, float& maxVDWrad) = 0;
   virtual std::list<AtomDescr> getAtDescOfAllPos(float &maxVDWrad) = 0;
   virtual void setHydAngle(double newAng, AtomPositions &xyz) = 0;

   const std::string& descr() const { return _descr; }
   void descr(const std::string& s) { _descr = s;    }

   double   bestScore() const { return _bestScore; }
   double   initScore() const { return _initScore; }

   bool bestHasBadBump() const { return _bestHasBadBump; }
   bool initHasBadBump() const { return _initHasBadBump; }

   int orientation(SearchStrategy ss=Mover::LOW_RES) const {
      return _orientation[ss!=Mover::LOW_RES];
   }
   int bestOrientation(SearchStrategy ss=Mover::LOW_RES) const {
      return _bestOrientation[ss!=Mover::LOW_RES];
   }

   void setBestOrientation(int io, double val, bool hasBadBump,
                           SearchStrategy ss=Mover::LOW_RES) {
      _bestOrientation[ss!=Mover::LOW_RES] = io;
      _bestScore      = val;
      _bestHasBadBump = hasBadBump;
   }

   virtual int  numOrientations(SearchStrategy ss=Mover::LOW_RES) const = 0;
   virtual void limitOrientations(bool, SearchStrategy ss=Mover::LOW_RES) = 0;
   virtual bool setOrientation(int oi, AtomPositions &xyz,
      SearchStrategy ss=Mover::LOW_RES) = 0;
   virtual bool isDefaultO(int oi, SearchStrategy ss=Mover::LOW_RES) const = 0;
   virtual std::string describeOrientation() const = 0;
   virtual double determineScore(AtomPositions &xyz,
      DotSphManager& dotBucket, int nBondCutoff,
      float probeRadius, float pmag, double& penalty,
      float &bumpScore, float &hbScore, bool& hasBadBump) = 0;

   bool valid() const { return _ok; }

   bool initOrientation(AtomPositions &xyz,
                           SearchStrategy ss=Mover::LOW_RES) {
      return setOrientation(0, xyz, ss);
   }
   bool nextOrientation(AtomPositions &xyz,
                           SearchStrategy ss=Mover::LOW_RES) {
      return (_orientation[ss!=Mover::LOW_RES]+1 < numOrientations(ss)) ?
	 setOrientation(_orientation[ss!=Mover::LOW_RES]+1, xyz, ss) : FALSE;
   }

   // connectivity within clique
   void addLink(const std::string &s, const MoverPtr lnk) {
	   _links.insert(std::make_pair(s, lnk));
   }
   bool linkExists(const std::string &s) const {
	   return _links.find(s) != _links.end();
   }
   int  numLinks() const { return _links.size(); }

   std::map<std::string, MoverPtr> links() const {
	   return _links;
   }

   double      flipMaxScore(int f) const { return   _flipMaxScore[f!=0]; }
   bool flipStateHasBadBump(int f) const { return _flipMaxBadBump[f!=0]; }

   void trackFlipMaxScore(int f, double val, bool hasBadBump);

   void makeNonAdjustable() { _ok = FALSE; _adj = FALSE; }
   
   virtual void dropBondedFromBumpingListForPDBrec( std::list< PDBrec * > & bumping, PDBrec* atom, int nBondCutoff ) const = 0;
protected:
   void initializeScoreIfNotSet(double val, bool hasBadBump) {
      if (!_initIsSet) {
	 _initScore      = val;
	 _initHasBadBump = hasBadBump;
	 _initIsSet      = TRUE;
      }
   }
   void rememberOrientation(int oi,
                           SearchStrategy ss=Mover::LOW_RES) {
      _orientation[ss!=Mover::LOW_RES] = oi;
   }
   void validateMemo()  { if (_adj) { _ok = TRUE; } }

   void resetFlipMaxScore(int f);
private:
	std::string      _descr;      // identifier

   double      _bestScore;
   double      _initScore;
   double      _flipMaxScore[2];   // max score for given state
   bool        _flipMaxBadBump[2]; // does the state have a bad bump?

   int         _orientation[2];
   int         _bestOrientation[2];

   bool        _bestHasBadBump; // does best orientation have a bad bump?
   bool        _initHasBadBump; // does best orientation have a bad bump?

   bool        _initIsSet;  // has the initial score been set?
   bool        _ok;         // is the memo valid?
   bool        _adj;        // is the memo adjustable?
   std::map<std::string, MoverPtr> _links; // connected items

   Mover(const Mover& m); // copy and assign not implemented
   Mover& operator=(const Mover& m);
};

void bondedList(const PDBrec& a, std::list<PDBrec*>& nearby, int nbnds,
		  std::list<PDBrec*>& atmList, std::list<PDBrec*>* bondedAtoms);
void countBonds(const PDBrec& src, const std::list<PDBrec*>& nearby,
	       int distcount, int maxcnt, std::list<PDBrec*>& atmList);
void resetMarks(std::list<PDBrec*>& lst);
bool visableAltConf(const PDBrec& a, bool onlyA);
bool interactingConfs(const PDBrec& a, const PDBrec& b, bool onlyA);
bool diffAltLoc(const PDBrec& a, const PDBrec& b);
int withinCovalentDist(const PDBrec& p, const PDBrec& q, double offset);
bool impossibleCovalent(const PDBrec& src, const PDBrec& targ, std::list<PDBrec*>& atmList);
double vdwGap(const PDBrec& p, const Point3d& pp,
              const PDBrec& q, const Point3d& qq);
bool foundInList(const PDBrec& a, const std::list<PDBrec*>& lst);
bool annularDots(const Point3d& dot, const PDBrec& src, const PDBrec& targ, float probeRadius);
double dot2srcCenter(const Point3d& dot, const PDBrec& src, const PDBrec& targ);
double kissEdge2bullsEye(float ra, float rb, float rp);

struct mpLess { 
	bool operator()(MoverPtr x, MoverPtr y) { 
		return x->descr() < y->descr();
	}
};

//bool mpLess (const MoverPtr& p1, const MoverPtr& p2);
//std::list<MoverPtr> mpSeqSort(const std::list<MoverPtr>& x);
//std::list<MoverPtr> mpSeqMerge(const std::list<MoverPtr>& x_list, const std::list<MoverPtr>& y_list);
//void mpSeqSplit(std::list<MoverPtr> x_list, std::list<MoverPtr>& y, std::list<MoverPtr>& z);
void mpSeqSort(std::list<MoverPtr>& x);
#endif
