// name: AtomPositions.h
// author: J. Michael Word
// date written: 8/1/97
// purpose: Interface for AtomPositions

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef ATOMPOSITIONS_H
#define ATOMPOSITIONS_H 1

#include <memory>
#include "PDBrec.h"
#include "ResBlk.h"
#include "DotSph.h"
#include "BumperPoint.h"
#include "Mover.h"
#include "CliqueList.h"
#include "neighbors.h"
#include "utility.h"
#include "AtomDescr.h"
#include "GraphToHoldScores.h"

extern bool UseSEGIDasChain; //jjh 130503, defined in reduce.cpp

class AtomPositions {
  private:
    class NullStream : public std::ostream {
        class NullBuffer : public std::streambuf {
        public:
            int overflow( int c ) { return c; }
        } m_nb;
    public:
        NullStream() : std::ostream( &m_nb ) {}
    };
    static NullStream nullStream;

  public:
     AtomPositions(int nblocks, bool onlyA, bool xplor, bool old, bool bbmodel, int nbCutoff,
        float minRegHBcut, float minChargedHBcut,
        float badBumpGapCut,
        DotSphManager& dotBucket, float probeRadius,
        float pmag, float occCutoff,
        bool verbose, bool showOrientScore,
        bool cliqueTicks, std::ostream& os = nullStream)
        : _onlyA(onlyA), _useXplorNames(xplor), _useOldNames(old), _bbModel(bbmodel),
    _nBondCutoff(nbCutoff),
    _min_regular_hb_cutoff(minRegHBcut),
    _min_charged_hb_cutoff(minChargedHBcut),
    _bad_bump_gap_cutoff(badBumpGapCut),
    _dotBucket(dotBucket),
    _probeRadius(probeRadius), _pmag(pmag),
    _occupancyCuttoff(occCutoff),
    _outputNotice(verbose),
    _showOrientScore(showOrientScore),
    _cliqueTicks(cliqueTicks), _os(os),
    scoreAtomsAndDotsInAtomsToScoreVector_( false ),
    scoreAtomsInAtomsInHighOrderOverlapList_( false ),
    atoms_to_score_ptr_(0), atoms_in_high_order_overlap_ptr_(0),
    _clqOfInt( false ) {}

  ~AtomPositions() {
  }

  int forceOrientations(const std::string& ofilename, std::vector<std::string>& notes);

   void put(std::shared_ptr<PDBrec> r) {
	   std::shared_ptr<PDBrec> temp = std::make_shared<PDBrec>(*r);
	   _xyzBlocks.insert(std::make_pair(LocBlk(r->loc()), temp));
   }

   // move atom to new xyz pos.
   void reposition(const Point3d& prev, const PDBrec& r);

   std::list< std::shared_ptr<PDBrec> > neighbors(const Point3d& p,
                         Coord mindist, Coord maxdist) const {
      return ::neighbors(p, mindist, maxdist, _xyzBlocks);
   }

   void insertRot(const PDBrec& hr, const PDBrec& c1,
                  const PDBrec& c2, const PDBrec& c3,
		  bool doOHSH, bool doNH3);

   void insertRotAromMethyl(const PDBrec& hr, const PDBrec& c1,
                  const PDBrec& c2, const PDBrec& c3); // for Arom methyls - Aram 08/13/12
   std::list<char> insertFlip(const ResBlk& rblk);
   void      insertFlip(std::shared_ptr<PDBrec> hr, std::list<char> alts_list);

   void doNotAdjust(const PDBrec& a);

   void finalizeMovers();

   CliqueList findCliques() const;

   int orientSingles(const std::list<MoverPtr>& singles);
   int orientClique(const std::list<MoverPtr>& clique, int limit); // returns -1 if abandoned
   int exhaustiveSearchOfClique(const std::list<MoverPtr>& clique);

   void describeChanges(std::list< std::shared_ptr<PDBrec> >& records,
	   std::list< std::shared_ptr<PDBrec> >::iterator& infoPtr, std::vector<std::string>& notes);

   int numChanges() const { return _motionDesc.size(); }

   void manageMetals(const ResBlk& rblk);

   void generateWaterPhantomHs(std::list< std::shared_ptr<PDBrec> >& waters);

   double atomScore(const PDBrec& a, const Point3d& p,
		float nearbyRadius, const std::list< std::shared_ptr<PDBrec> >& exclude,  // JSS: no need to copy
		const DotSph& dots, float prRadius, bool onlyBumps,
		float &bumpSubScore, float &hbSubScore, bool &hasBadBump);

	//double determineScoreForMoverIn3WayOverlap(
	//	std::list< AtomDescr > * atomsIn3WayOverlap,
	//	MoverPtr mover
	//);

	void CollectBumping(const AtomDescr& ad, std::list< std::shared_ptr<PDBrec> >& bumping);
   float & getMaxFoundVDWRad() { return _maxVDWFound;}
   int getNBondCutoff() const {return _nBondCutoff;}

   float determineScoreForMover(
   	MoverPtr mover,
   	std::vector< std::pair< AtomDescr, DotsForAtom * > >  & atoms_to_score,
   	double & penalty
   );

   //float scoreMoverInHighOrderOverlap(
	//	std::list< AtomDescr > & atomsInHighOrderOverlap,
	//	MoverPtr mover );

	bool outputNotice() const {return _outputNotice;}

private:
   //AtomPositions(const AtomPositions& a);            // can't copy
   //AtomPositions& operator=(const AtomPositions& a); // can't assign

   std::multimap<LocBlk, std::shared_ptr<PDBrec> > _xyzBlocks;
   std::map<std::string,  MoverPtr>       _motionDesc;
   std::list< std::shared_ptr<BumperPoint> >          _excludePoints;

   const bool                _onlyA;
   const bool                _useXplorNames;
   const bool                _useOldNames;
   const bool		     _bbModel;
   const int                 _nBondCutoff;
   const float               _min_regular_hb_cutoff;
   const float               _min_charged_hb_cutoff;
   const float               _bad_bump_gap_cutoff;
   DotSphManager&            _dotBucket;
   const float               _probeRadius;
   const float               _pmag;
   const float               _occupancyCuttoff;
   const bool                _outputNotice;
   const bool                _showOrientScore;
   const bool                _cliqueTicks;
   std::ostream&             _os;
	float	_maxVDWFound;

   bool
	initializeCliqueMovers(
		std::list< MoverPtr > const & clique,
		std::vector< MoverPtr > & item,
		int const numItems
	);

	void
	setNumStatesForNodes(
		std::vector< MoverPtr > & item,
		int const numItems,
		std::vector< int > & num_states,
		std::vector< std::vector< float > > & penalties
	);


private:

	bool scoreAtomsAndDotsInAtomsToScoreVector_;
	bool scoreAtomsInAtomsInHighOrderOverlapList_;
	std::vector< std::pair< AtomDescr, DotsForAtom * > > * atoms_to_score_ptr_;
	std::list< AtomDescr > * atoms_in_high_order_overlap_ptr_;

	int  _clqDots;				//ANDREW: for keeping track of how many dots go into the score for a network
	int  _bgDots;				//ANDREW: used only when the flag DEBUGDOTCOUNTS is set.
	bool _clqOfInt;			//ANDREW: this will let me detect once whether the clique I'm examining is the one I'm interested in
									//and keep that fact around when the scoring is being done.


public:
   int SearchClique(std::list<MoverPtr> clique, int time_limit);

};
#endif
