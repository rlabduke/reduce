// Name: reduce.h
// Author: J. Michael Word
// Date Written: 7/15/97
// Purpose: add hydrogens to a Protein Data Bank file and
//          writing a new Protein Data Bank file to standard output
//
// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999-2016 J. Michael Word
// **************************************************************

#include <iostream>
#include <list>
#include <vector>
#include <memory>
#include "PDBrec.h"
#include "CTab.h"
#include "AtomPositions.h"

struct SummaryStats {
  SummaryStats() : _H_found(0), _H_HET_found(0),
    _H_removed(0), _H_HET_removed(0),
    _H_added(0), _H_HET_added(0),
    _H_standardized(0), _H_HET_standardized(0),
    _num_atoms(0), _conect(0),
    _num_adj(0), _num_renamed(0) {}

  int _H_found, _H_HET_found;
  int _H_removed, _H_HET_removed;
  int _H_added, _H_HET_added;
  int _H_standardized, _H_HET_standardized;
  int _num_atoms, _conect;
  int _num_adj, _num_renamed;
};
extern SummaryStats Tally;

/// @brief Read the PDB records from the specified input stream.
/// @param [in] s String with one PDB record per line.
/// @return List (one entry per model) of PDB records found in the file.
///       Each list contains all of the records that are outside of
///       MODEL...ENDMDL records along with the records for the model
///       that was read, with the first one being MODEL 1.
extern std::vector< std::list< std::shared_ptr<PDBrec> > > inputModels(std::string s);

/// @brief Drop hydrogens from a model in place.
/// @param [in] RemoveATOMHydrogens Should we remove hydrogens in ATOM records?
/// @param [in] RemoveOtherHydrogens Should we remove hydrogens in non_ATOM records?
extern void dropHydrogens(std::list< std::shared_ptr<PDBrec> >& records,
  bool RemoveATOMHydrogens, bool RemoveOtherHydrogens);

/// @brief Optimize all of the records passed in in place.
/// @return 0 on success, 1 on abandoned due to too many permutations.
extern int optimize(AtomPositions &xyz, std::vector<std::string> &adjNotes);

//SJ 08/03/2015 for printing all models together
extern void outputRecords_all(
  const std::vector<std::list< std::shared_ptr<PDBrec> > >& all_records,
  std::ostream& os = std::cout);

/// @brief Output all records into a string with one record per line.
///
/// This reassembles the multiple models into a single combined file that
/// has the boiler-plate information once and all models in the same file.
/// It also reports the version of version of Reduce that was run in a
/// USER MOD line.
/// @param all_records vector of models to combine.
/// @return String with one record described on each line
extern std::string outputRecords_all_string(
    const std::vector<std::list< std::shared_ptr<PDBrec> > >& all_records);

/// @brief Check the list of PDB records to see if we should use segment ID as chain
/// @return TRUE if we should use the segment ID as the chain, FALSE if not
extern bool checkSEGIDs(std::list< std::shared_ptr<PDBrec> >& rlst);

extern void scanAndGroupRecords(std::list< std::shared_ptr<PDBrec> >& rlst, AtomPositions& xyz,
  std::list< std::shared_ptr<PDBrec> >::iterator& startAtoms);
extern void reduceList(CTab& db, std::list< std::shared_ptr<PDBrec> >& records,
  AtomPositions& xyz, std::vector<std::string>& fixNotes);

/// Library-scoped global variables that are in the global name space.
/// @todo Move these out of the global name space
extern const char *versionString;
extern const char *shortVersion;
extern const char *referenceString;
extern const char *electronicReference;

extern bool Verbose;    // do we write processing notes to stdout?
extern bool KeepConnections;
extern bool StandardizeRHBondLengths;
extern bool ProcessConnHydOnHets;
extern bool BuildHisHydrogens;
extern bool SaveOHetcHydrogens;
extern bool UseXplorNames;
extern bool UseOldNames;
extern bool BackBoneModel;
extern bool DemandRotAllMethyls;
extern bool RotExistingOH;
extern bool NeutralTermini;
extern bool DemandRotNH3;
extern bool DemandRotExisting;
extern bool DemandFlipAllHNQs;
extern bool DoOnlyAltA; //jjh changed default 111118
extern bool OKProcessMetMe; //cjw changed default 160602
extern bool OKtoAdjust;
extern bool ShowCliqueTicks;
extern bool ShowOrientScore;
extern bool StringInput;
extern bool ShowCharges;
extern bool UseNuclearDistances; //jjh 130111
extern bool UseSEGIDasChain; //jjh 130503
extern bool ProcessedFirstModel; //jjh 130718
extern bool RenameFlip; // SJ - 09/25/2015 - flag to specify if the final PDB file will have the coordinates according to
                        //the rename atoms flip (is the flag is TRUE) or the new rot hinge dock flip (if the flag is FALSE, this is default). Can be set to true by the -renameflip flag in the commandline.
extern bool GenerateFinalFlip; // SJ - 09/04/2015 to keep track of when scoring and decision of flips finishes and when the final PDB
                                       //coordinates are being generated. This is set to true after all the calculations are done, unless the RenameFlip flag is TRUE

extern int MaxAromRingDih;   // max dihedral angle in planarity check for aromatic rings  120724 - Aram

extern int MinNTermResNo;   // how high can a resno be for n-term?
extern int NBondCutoff;   // how many bonds away do we drop?
extern int ExhaustiveLimit;  //time limit, in seconds, to spend in brute force enumeration for a single clique
extern float ProbeRadius; // how big is the probe in VDW calculations?
extern float VdwDotDensity; // how many dots per sq Angstroms in VDW calculations?
extern float OccupancyCutoff;// lowest occupancy considered when determining score
extern float WaterBcutoff; // limit for water B values
extern float WaterOCCcutoff;// limit for water occupancy
extern float PenaltyMagnitude;// score bias towards original orientation (changed from 0.0 in 2.13.0)
extern float MinRegHBgap; // Hbonds with greater gaps start to bump
extern float MinChargedHBgap; // charged Hbonds start to bump at this point
extern float BadBumpGapCut; // bump is bad if >= than this
extern float NonMetalBumpBias;//bumps if H closer than atom radius, plus this
extern float MetalBumpBias;// ditto, for metals
extern float GapWidth; // half width for detecting chain breaks between residues
                       // (center at 1.4; default allow 1.1-1.7 for accepting connected residues)

extern std::string OFile; // if file exists, given orientations forced
extern bool UseSEGIDtoChainMap; // if true, override some chain ids
extern bool StopBeforeOptimizing;  // If true, handle hydrogen drops/adds and then quit
extern bool AddWaterHydrogens;	  // If true, add phantom hydrogens to waters
extern bool AddOtherHydrogens;	  // If true, add hydrogens to non-water atoms
extern bool RemoveATOMHydrogens; // If true, remove hydrogens from ATOM records
extern bool RemoveOtherHydrogens;// If true, remove hydrogens from non-ATOM records
