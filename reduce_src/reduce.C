// Name: reduce.C
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
// Copyright (C) 1999-2003 J. Michael Word
// **************************************************************
//
// History: this new version of reduce has been completely re-written
//          to include het groups and rotations
//
// 10/ 1/97 - jmw - changed rotation score to dot score, modified geometry,
//                  of CH2 and mcNH angles to be compatible with ECEPP,
//                  and added option to permit all methyls to rotate.
// 10/22/97 - jmw - fixed bug in rotation code which did not separate
//                  alternate conformations
// 11/12/97 - jmw - support for flips of HIS and ASN/GLN residues,
//                  orientable waters and Hbonds to aromatics
//                  (includes some analysis of pairing of flips and rots)
// 12/10/97 - jmw - fixed waters to have phantom H atoms,
//                  better rotational searching, metal covalent radii
//  2/ 7/98 - jmw - broke out motions and re-wrote search,
//                  use only A conf, fix penalty
//  2/18/98 - jmw - fixup ambiguous atom names, no creation of altB Hs,
//                  MTO waters, added -BUILD and other flags,
//                  better recognise S-S or C-O-C bonds, etc.
//  3/ 1/98 - jmw - fix orientations
//  3/25/98 - jmw - updated metal radii, added HET dict env var,
//                  expanded fixed orientations, changed -OH default
//  4/ 8/98 - jmw - fixed bug, occupancy of water phantomHs now > 0.0 !
//  4/18/98 - jmw - bad bump check, short water h, new limit on HB, ...
//  4/26/98 - jmw - new hires search strategy
//  4/29/98 - jmw - fixed min rotation angle and penalty bugs
//  5/13/98 - jmw - fixed another min rotation angle bug where the coarse score
//                  is the best that can be obtained. Also updated output to
//                  the header to document K/C/X/F categories.
//  5/15/98 - jmw - fine tuned cases where group is fixed by metal or modification
//  7/ 3/98 - jmw - fixed bug: ASN/GLN sc C=O carbon radii was 1.75, now 1.65
//  8/ 7/98 - jmw - stop putting notes in the segID, elem & charge fields
//                  because PDB use of these fields are now being expanded
//                  by naming convention differences with XPLOR and XPLORs
//                  use of SEGID rather than CHAINID
//  8/12/98 - jmw - changed how we figure num of permutations
//  8/14/98 - jmw - now search for het dictionary in current directory
//  9/ 1/98 - jmw - worked on portability by cleaning up g++ warnings
//  9/13/98 - jmw - figured out a work-around to g++ template and
//                  static const class initializer problems
//  1/ 7/99 - jmw - consolidate Linux and Mac changes with SGI src
//  3/16/99 - jmw - extended atom name parsing for wierd HETs with col1 ABCDEFGs
//  4/ 6/99 - jmw - added OW to list of names for oxygen in water
//  7/ 7/99 - jmw - Improved portability (near->nearr, List <T>::, Makefiles),
//                  made penalty 1.0, added reference, changed basic io hooks
//  9/ 2/99 - jmw - Updated Makefile list and Utility.h for DEC alpha support
// 10/ 5/99 - jmw - Updated a Makefiles and the main in reduce.C for sgi6.5 compiler
//  8/ 4/00 - jmw - Added -segid flag to support segment identifiers
// 10/30/00 - jmw - Modified Seq mergesort const/non const stuff
//  4/19/01 - jmw - Added support for left justified A/T/U/C/G for nucleic acids
//                  which fixed a rare but nasty bug in the water recognition
//  5/ 7/01 - jmw - Stopped output with -quiet while fixing orientation
//                  and finally dealt with string literal conversion messages
//  5/31/01 - jmw - Pass signal of abandoned clique search in rc,
//                  added new flag to allow mapping of segids to chains,
//                  changed properties of Nterminal fragment nitrogens to acceptor
// 10/ 4/01 - jmw - fiddled to get compiled on RH linux 7 (gcc 2.96)
//  5/24/02 - jmw - added control over hbump
//  1/23/03 - jmw -v2.16 - Changed global MinChargedHBgap: 0.4 => 0.8
//                         Changed phosphorus properties to drop acceptor status
//  3/06/03 - jmw -v2.17 - Cleaning declarations found by Leo and Andrew such as
//                         adding compile flag -DOLD_CXX_DEFNS to select const longs
//                         instead of std::_Ios_Fmtflags in AtomPositions.C
//  3/31/03 - jmw -v2.18 - Fixed spelling mistake by renaming cuttoff to cutoff throughout.
//                         Edited help for -H2OBcutoff# to suggest an integer value
//                         Made -H2OOCCcutoff#.# parse a real (was integer)
//  4/ 3/03 - jmw -v2.19 - Moved -Help to the end of the command line parsing.
//  6/ 3/03 - jmw -v2.20 - updated to isoC++ style std includes: #include <cstring>,
//                         changed String class to Stringclass class,
//                         nice speedups from Jack Snoeyink
//  6/ 4/03 - jmw -v2.21 - added out for OLD_STD_HDRS for sgi plus other sgi polishing
// 11/ 7/03 - jmw -v2.22 - fixed bug with creating hydrogens at the end of a triple bond
//                         and supressed the warnings about connecting nucleic acid bases
//
// 040509dcr reduce.C ln 455: NO renumber, 1459: NO RXR msg
//
// 041113dcr reduce.C reconstructed main to loop over any NMR models in the file
//           a specific model can still be specified.
//
//changes
//
// 050314 dcr incorporated Mike's version of 030703, to wit:
//
//  6/15/06 - apl -v3.0  - incorporated decomposition of scoring function into interactions of
//                         atom singles, atom pairs, atom tripples and all other
//                         higher-order atom interactions.  incorporated
//                         dynamic programming. incorporated changes from v2.21.mod_dcr
//                         disabling hydrogen atom numbering and skipInfo. removing
//                         PDBrecNAMEout as in v2.999. Disabling hydrogen bonds to his
//                         aromatic carbon atoms.
//  6/19/06 - apl -v3.01- decomposing the scoring function in terms of which dots should
//                        be scored with which hyperedges.  Incorporating S3 reduction rules
//                        into dynamic programming.  Incorporating code to handle 4-way overlap
//                        (but not five way overlap) though I have not observed any 4-way
//                        overlap using the new decomposition scheme (it would show up in the previous scheme).
//  6/24/06 - apl -v3.02- incorporating the additions to v2.23 that deal with multiple NMR models
//                        in a single file.  Altering output from dp to be less intrusive and a
//                        little more informative.
//                        Adding new reduction rule for vertices with exactly 1 state (i.e. no
//                        real options) which speeds up dynamic programming for the second-round
//                        of optimizations that calculate the optimal network states for sub-optimal
//                        flip states.  GLN and ASN will have 1 state in these optimizations.
// 7/03/06 - apl -      - Fixing USER MOD Set ordering in output PDB.  Fixing "flip" records in columns
//                        82 to 85 for the atoms on flipped residues.  Fixing atom placement plan bug
//                        in genHydrogens() that manifested itself as aberant behavior when presented
//                        with several (3 or more) alternate conformations.
// 7/09/06 - apl -      - Adding #include <cassert> for gcc3.3.4 builds
// 7/11/06 - apl -      - Fixing "node with 0 states" bug.

#pragma warning(disable:4786) 
#pragma warning(disable:4305) 
#pragma warning(disable:4800) 

static const char *versionString =
     "reduce: version 3.02  7/09/06, Copyright 1997-2006, J. Michael Word";

static const char *shortVersion    = "reduce.3.02.060709";
static const char *referenceString =
                       "Word, et. al. (1999) J. Mol. Biol. 285, 1735-1747.";
static const char *electronicReference = "http://kinemage.biochem.duke.edu";

#ifdef OLD_STD_HDRS
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <string.h>
#include <ctype.h>
#else
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cctype>
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#endif

#ifndef NOSYSSTATS
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif
#include <vector>
#include <list>
#include "utility.h"
#include "CTab.h"
#include "StdResH.h"
#include "ResBlk.h"
#include "pdb++.h"
#include "PDBrec.h"
#include "AtomPositions.h"

#define DIRECTORY_SEP_CHAR '/'

int ReturnCodeGlobal = 0;
#define ABANDONED_RC 1

bool Verbose = TRUE;    // do we write processing notes to stdout?
bool KeepConnections          = TRUE;
bool StandardizeRHBondLengths = TRUE;
bool ProcessConnHydOnHets     = TRUE;
bool RemoveHydrogens          = FALSE;
bool BuildHisHydrogens        = FALSE;
bool SaveOHetcHydrogens       = TRUE;
bool UseXplorNames            = FALSE;
bool DemandRotAllMethyls      = FALSE;
bool RotExistingOH            = FALSE;
bool DemandRotNH3             = TRUE;
bool DemandRotExisting        = FALSE;
bool DemandFlipAllHNQs        = FALSE;
bool DoOnlyAltA               = TRUE;
bool OKProcessMetMe           = TRUE;
bool OKtoAdjust               = TRUE;
bool ShowCliqueTicks          = TRUE;
bool ShowOrientScore          = FALSE;

int MinNTermResNo     = 1;   // how high can a resno be for n-term?
int ModelToProcess    = 1;   // which model to work on, 
                             // >0 is a model to work on  041113
int ModelSpecified    = 0;   // commandline model specified  041113
int ModelNext         = 0;   // next model to process  041113
int ModelActive       = 0;   // found the next model and working on it  041113
int NBondCutoff       = 3;   // how many bonds away do we drop?
int ExhaustiveLimit   = 600;  //time limit, in seconds, to spend in brute force enumeration for a single clique
float ProbeRadius     = 0.0; // how big is the probe in VDW calculations?
float VdwDotDensity   =16.0; // how many dots per sq Angstroms in VDW calculations?
float OccupancyCutoff = 0.01;// lowest occupancy considered when determining score
float WaterBcutoff    =40.0; // limit for water B values
float WaterOCCcutoff  = 0.66;// limit for water occupancy
float PenaltyMagnitude= 1.00;// score bias towards original orientation (changed from 0.0 in 2.13.0)
float MinRegHBgap     = 0.6; // Hbonds with greater gaps start to bump
float MinChargedHBgap = 0.8; // charged Hbonds start to bump at this point
float BadBumpGapCut   = 0.4; // bump is bad if >= than this
float NonMetalBumpBias= 0.125;//bumps if H closer than atom radius, plus this
float MetalBumpBias   = 0.865;// ditto, for metals

std::string OFile; // if file exists, given orientations forced
bool UseSEGIDtoChainMap = FALSE; // if true, override some chain ids

#ifndef HET_DICTIONARY
#define HET_DICTIONARY "reduce_het_dict.txt"
#endif
std::string DBfilename( HET_DICTIONARY );

enum ConnType {NTERM_RES, CONNECTED_RES, FRAGMENT_RES};

struct SummaryStats {
   SummaryStats() : _H_found(0),        _H_HET_found(0),
                    _H_removed(0),      _H_HET_removed(0),
                    _H_added(0),        _H_HET_added(0),
                    _H_standardized(0), _H_HET_standardized(0),
		    _num_atoms(0),      _conect(0),
		    _num_adj(0),        _num_renamed(0) {}

   int _H_found,        _H_HET_found;
   int _H_removed,      _H_HET_removed;
   int _H_added,        _H_HET_added;
   int _H_standardized, _H_HET_standardized;
   int _num_atoms,      _conect;
   int _num_adj,	_num_renamed;
};

SummaryStats Tally;

char* parseCommandLine(int argc, char **argv);
void processPDBfile(istream& ifs, char *pdbFile, ostream& ofs);
void establishHetDictionaryFileName(void);
void reduceHelp(bool showAll);
istream& inputRecords(istream& is, std::list<PDBrec*>& records);
ostream& outputRecords(ostream& os, const std::list<PDBrec*>& records);
void dropHydrogens(std::list<PDBrec*>& records);
void reduceList(CTab& db, std::list<PDBrec*>& records,
				AtomPositions& xyz, std::vector<std::string>& fixNotes);
void scanAndGroupRecords(std::list<PDBrec*>& rlst, AtomPositions& xyz,
						 std::list<PDBrec*>::iterator& startAtoms);
void renumberAndReconnect(std::list<PDBrec*>& rlst);
void renumberAtoms(std::list<PDBrec*>& rlst);
void analyzeRes(CTab& db, ResBlk* pb, ResBlk* cb, ResBlk* nb,
				AtomPositions& xyz, std::list<PDBrec*>& waters, 
				std::vector<std::string>& fixNotes, std::list<PDBrec*>& rlst);
bool isConnected(ResBlk* a, ResBlk* b);
void genHydrogens(const atomPlacementPlan& pp, ResBlk& theRes, bool o2prime,
				  AtomPositions& xyz, std::list<char>& resAlts,
		  std::vector<std::string>& fixNotes, std::list<PDBrec*>& rlst);
void noteWaterInfo(ResBlk& r, std::list<PDBrec*>& waters);
void findAndStandardize(const char* name, const char* resname, ResBlk& theRes);
void stdBondLen(float dist, PDBrec& ourHydrogen, std::list<PDBrec*>& firstAtoms,
                                     const ElementInfo& e);
void fixupHeavyAtomElementType(ResBlk& theRes, CTab& hetDB);
void noteAddedInTally(const PDBrec& theH);
void noteStdInTally(const PDBrec& theH);
void noteRemovedInTally(const PDBrec& theH);
int findConnAtoms(const atomPlacementPlan& pp, ResBlk& theRes, char hac,
   PDBrec& r0atom, PDBrec& r1atom, PDBrec& r2atom);
bool okToPlaceHydHere(const PDBrec& theHatom,
   const atomPlacementPlan& pp, const PDBrec& a1, const PDBrec& a2,
   AtomPositions& xyz, bool& doNotAdjustSC,
   std::vector<std::string>& fixNotes);
void recordSkipInfo(bool skipH, std::vector<std::string>& fixNotes,
   const PDBrec& theHatom, const PDBrec& heavyAtom,
   std::list<PDBrec*>& nearr, const char * msg);

//#ifdef MACAPPL
//   char *getOptionsOnTheMac();

//int main() {
//   Verbose = FALSE;
//   ShowCliqueTicks = FALSE;
//
//   char *pdbFile = getOptionsOnTheMac();
//
//   ifstream theinputstream(pdbFile);
//   processPDBfile(theinputstream, pdbFile, ofstream("dump.out"));
//
//   return ReturnCodeGlobal; // one pass and then we quit
//}
//#else

int main(int argc, char **argv) {

   char *pdbFile = parseCommandLine(argc, argv);

   if(pdbFile)
   {/*can do multiple passes by rewinding input file*/
      //ifstream theinputstream(pdbFile); //would need rewind
      while(ModelToProcess) /* 041113 */
      {
         ifstream theinputstream(pdbFile); //declare each time, avoid rewind
         processPDBfile(theinputstream, pdbFile, cout); 
         if(ModelSpecified) {ModelToProcess = 0;} /* did it, so quit */
         else if(ModelNext > 0) 
         { 
            ModelToProcess = ModelNext;
            ModelNext = 0; /*perhaps to be rediscovered in PDB file*/
            //theinputstream::rewind; /*theinputstream undeclared*/
//cerr<<"about to rewind"<<endl;
//            rewind(theinputstream); /*gives warnings*/
//cerr<<"just did rewind"<<endl;
         }
         else {ModelToProcess = 0;} /*ModelNext==0, time to quit*/
      }
   }
   else
   {/*presume stdin for pdb file info*/
      processPDBfile(cin,            pdbFile, cout); 
   }
   return ReturnCodeGlobal;
}
//#endif

void processPDBfile(istream& ifs, char *pdbFile, ostream& ofs) {
   if (Verbose) {
      cerr << versionString << endl;
      if (pdbFile) {
	 cerr << "Processing file: \"" << pdbFile << "\"" << endl;
      }
      else {
	 cerr << "Processing file: --standard input--" << endl;
      }
      if (ProcessConnHydOnHets && ! RemoveHydrogens) {
	 cerr << "Database of HETATM connections: \"" << DBfilename << "\"" << endl;
      }
      if (ModelToProcess != 1) {
	 cerr << "Using Model: \"" << ModelToProcess << "\"" << endl;
      }
   }

   std::list<PDBrec*> records;

   inputRecords(ifs, records);       // read all the PDB records in the file

   if (RemoveHydrogens) {
      dropHydrogens(records);

      if (Verbose) {
	 cerr << "Trimming: removed " <<Tally._H_removed << " hydrogens ("
	    << Tally._H_HET_removed << " hets)" << endl;
      }

      ofs << "USER  MOD "<< shortVersion <<" removed " << Tally._H_removed
           << " hydrogens (" << Tally._H_HET_removed << " hets)" << endl;
   }
   else {
	 if (Verbose) {
	    if (! OFile.empty()) {
	       cerr << "Using orientation info in \""<<OFile<<"\"." << endl;
	    }
	    if (UseSEGIDtoChainMap) {
               PDBrec::DumpSEGIDtoChainMap(cerr, "Mapping ");
	    }

	    if (DoOnlyAltA) {
	       cerr << "Processing only 'A' conformations." << endl;
	    }
	    cerr << "VDW dot density = " << VdwDotDensity << "/A^2" << endl;
	    if (ProbeRadius > 0.01) {
	       cerr << "Probe radius = " << ProbeRadius << "A" << endl;
	    }
	    if (PenaltyMagnitude > 0.0001 || PenaltyMagnitude < -0.0001) {
	       cerr << "Orientation penalty scale = "<< PenaltyMagnitude
	            <<" (" << int(PenaltyMagnitude*100.0+0.5) << "%)" << endl;
	    }
	    else {
	       cerr << "No Orientation penalty (-penalty0.0)" << endl;
	    }
	    cerr << "Eliminate contacts within " << NBondCutoff << " bonds." << endl;
	    cerr << "Ignore atoms with |occupancy| <= "
	         << OccupancyCutoff << " during adjustments." << endl;
	    cerr << "Waters ignored if B-Factor >= " << WaterBcutoff
	                    << " or |occupancy| < " << WaterOCCcutoff<< endl;
#ifdef AROMATICS_ACCEPT_HBONDS
	    cerr << "Aromatic rings in amino acids accept hydrogen bonds." << endl;
#endif
	    if (BuildHisHydrogens) {
	       cerr << "Building His ring NH Hydrogens." << endl;
	    }
	    if (DemandFlipAllHNQs) {
	       cerr << "Flipping Asn, Gln and His groups." << endl;
	       cerr << "For each flip state, bumps where gap is more than "
	         << BadBumpGapCut << "A are indicated with \'!\'." << endl;

	    }
	    if (SaveOHetcHydrogens) {
	       cerr << "Building or keeping OH & SH Hydrogens." << endl;
	       if (RotExistingOH && !DemandRotExisting) {
		  cerr << "Rotating existing OH & SH Hydrogens" << endl;
	       }
	    }
	    else {
	       cerr << "Dropping existing OH & SH Hydrogens." << endl;
	    }
	    if (DemandRotNH3) {
	       cerr << "Rotating NH3 Hydrogens." << endl;
	    }
	    if (DemandRotAllMethyls) {
	       cerr << "Rotating ALL methyls." << endl;
	    }
	    if (OKProcessMetMe) {
	       if (! DemandRotAllMethyls) {
		  cerr << "Processing Met methyls." << endl;
	       }
	    }
	    else {
	       cerr << "Not processing Met methyls." << endl;
	    }
	    if (DemandRotExisting) {
	       cerr << "Rotating pre-existing groups." << endl;
	    }
	 }

      CTab hetdatabase(DBfilename, 1000);

      DotSphManager dotBucket(VdwDotDensity);

      AtomPositions xyz(2000, DoOnlyAltA, UseXplorNames, NBondCutoff,
                        MinRegHBgap, MinChargedHBgap, BadBumpGapCut,
			dotBucket, ProbeRadius,
			PenaltyMagnitude, OccupancyCutoff,
			Verbose, ShowOrientScore, ShowCliqueTicks, cerr);

//      NonConstListIter<PDBrec> infoPtr(records); // info on changes can be
                                                 // inserted before infoPtr
	  std::list<PDBrec*>::iterator infoPtr = records.begin();

      scanAndGroupRecords(records, xyz, infoPtr);

// if some sidechain needs adjustment...
	  std::vector<std::string> adjNotes;
      Tally._num_adj = 0;

      reduceList(hetdatabase, records, xyz, adjNotes);

      if ((OKtoAdjust || ! OFile.empty()) && xyz.numChanges() > 0) {
	 xyz.finalizeMovers();
      }
      if (! OFile.empty()) {
	 Tally._num_adj += xyz.forceOrientations(OFile, adjNotes);
      }

      if (OKtoAdjust && xyz.numChanges() > 0) {

	 CliqueList clst = xyz.findCliques();

	 if (Verbose) { clst.describe(cerr); }

// adjust singletons

	 Tally._num_adj += xyz.orientSingles(clst.singles());

// adjust cliques

	 std::list< std::list<MoverPtr> > cc_list = clst.cliques();
	//cerr << "start: " << cc_list.size() << endl;
	 for (std::list< std::list<MoverPtr> >::iterator cc = cc_list.begin(); cc != cc_list.end(); ++cc) {
		//cerr << "start2" << endl;
	    int nscnt = xyz.orientClique(*cc, ExhaustiveLimit);
		if (nscnt > 0) { Tally._num_adj += nscnt; }
	    else { // too many permutations, make note
	      ReturnCodeGlobal = ABANDONED_RC;
	    }
	 }

// record clique and singleton adjustments

	 for (int jj = 0; jj < clst.numCliques(); jj++) {
	    clst.formatClique(adjNotes, jj);
	 }
	 clst.formatSingles(adjNotes);
	 xyz.describeChanges(records, infoPtr, adjNotes);
      }

      if (Verbose) {
	 cerr << "Found " <<Tally._H_found << " hydrogens ("
	      << Tally._H_HET_found << " hets)" << endl;
	 cerr << "Standardized " <<Tally._H_standardized << " hydrogens ("
	      << Tally._H_HET_standardized << " hets)" << endl;
	 cerr << "Added " <<Tally._H_added << " hydrogens ("
	      << Tally._H_HET_added << " hets)" << endl;
	 cerr << "Removed " <<Tally._H_removed << " hydrogens ("
	      << Tally._H_HET_removed << " hets)" << endl;
	 if (Tally._num_adj > 0) {
	    cerr << "Adjusted " << Tally._num_adj << " group(s)" << endl;
	 }
	 if (Tally._num_renamed > 0) {
	    cerr << "Renamed and marked " << Tally._num_renamed
	         << " ambiguous 'A' atom name(s)" << endl;
	 }
      }

      ofs << "USER  MOD "<< shortVersion
           <<" H: found="<<Tally._H_found
           <<", std="<< Tally._H_standardized
           << ", add=" << Tally._H_added
           << ", rem=" << Tally._H_removed
           << ", adj=" << Tally._num_adj << endl;
      if (Tally._num_renamed > 0) {
	 ofs << "USER  MOD renamed " << Tally._num_renamed
	      << " ambiguous 'A' atoms (marked in segID field)"<< endl;
      }
      if (UseSEGIDtoChainMap) {
         PDBrec::DumpSEGIDtoChainMap(ofs, "USER  MOD mapped ");
      }
   }

   if (Tally._H_added || Tally._H_removed) {
      ;//renumberAndReconnect(records);
   }

   // write out the results...
   outputRecords(ofs, records);

   if (Verbose) {
      cerr << "If you publish work which uses reduce, please cite:"
           << endl << referenceString << endl;
      cerr << "For more information see " << electronicReference << endl;
   }

   std::for_each(records.begin(), records.end(), DeleteObject());
}

void establishHetDictionaryFileName(void) {
	int i =DBfilename.find_last_of(DIRECTORY_SEP_CHAR); 
	std::string localfile = DBfilename.substr(i+1, DBfilename.length());
#ifdef NOSYSSTATS
	const int rc = 1; // force the name on the MAC
#else
	struct stat filestatbuf;
	const int rc = stat(localfile.c_str(), &filestatbuf);
#endif
	if (rc == 0) {
		DBfilename = localfile;
	}
	else {
		const char *hetdbAlias = getenv("REDUCE_HET_DICT");
		if (hetdbAlias && (strlen(hetdbAlias) > 0)) {
			DBfilename = hetdbAlias;
		}
	}
}

#ifdef MACAPPL
char *getOptionsOnTheMac() {
   establishHetDictionaryFileName(); // this may need modifying

   BuildHisHydrogens  = TRUE; // this block is the same as the old -build
   SaveOHetcHydrogens = TRUE;
   RotExistingOH      = TRUE;
   DemandFlipAllHNQs  = TRUE;

   return "test.pdb";
}
#endif

char* parseCommandLine(int argc, char **argv) {
   char *pdbfile = NULL;
   int nfile = 0, n;

   establishHetDictionaryFileName();

   if (argc <= 1) {
      reduceHelp(FALSE);
   }
   for (int i = 1; i < argc; i++) {
      char *p = argv[i];
      if (p[0] == '-') {
	 if (p[1] == '\0') {
	    nfile = 1;
	    pdbfile = NULL; // i.e. standard input
	 }
	 else if(compArgStr(p+1, "BUILD", 5)){
	    BuildHisHydrogens  = TRUE;
	    SaveOHetcHydrogens = TRUE;
	    RotExistingOH      = TRUE;
	    DemandFlipAllHNQs  = TRUE;
	 }
	 else if(n = compArgStr(p+1, "Quiet", 1)){
	    Verbose = FALSE;
	 }
	 else if(n = compArgStr(p+1, "NOTICKs", 6)){
	    ShowCliqueTicks = FALSE;
	 }
	 else if(n = compArgStr(p+1, "SHOWSCore", 6)){
	    ShowOrientScore = TRUE;
	 }
	 else if(n = compArgStr(p+1, "NOCon", 3)){
	    KeepConnections = FALSE;
	 }
	 else if(n = compArgStr(p+1, "NOROTMET", 8)){
	    OKProcessMetMe = FALSE;
	 }
	 else if(n = compArgStr(p+1, "NOADJust", 5)){
	    OKtoAdjust = FALSE;
	 }
	 else if(n = compArgStr(p+1, "HIS", 3)){
	    BuildHisHydrogens = TRUE;
	 }
	 else if(n = compArgStr(p+1, "OH", 2)){
	    SaveOHetcHydrogens = TRUE;
	 }
	 else if(n = compArgStr(p+1, "NOOH", 4)){
	    SaveOHetcHydrogens = FALSE;
	 }
	 else if(n = compArgStr(p+1, "Xplor", 1)){
	    UseXplorNames = TRUE;
	 }
	 else if(n = compArgStr(p+1, "Trim", 1)){
	    RemoveHydrogens = TRUE;
	 }
	 else if(n = compArgStr(p+1, "Keep", 1)){
	    StandardizeRHBondLengths = FALSE;
	 }
	 else if(n = compArgStr(p+1, "ALLMEthyls", 5)){
	    DemandRotAllMethyls = TRUE;
	 }
	 else if(n = compArgStr(p+1, "ROTEXist", 5)){
	    DemandRotExisting = TRUE;
	 }
	 else if(n = compArgStr(p+1, "ROTNH3", 6)){
	    DemandRotNH3 = TRUE;
	 }
	 else if(n = compArgStr(p+1, "NOROTNH3", 8)){
	    DemandRotNH3 = FALSE;
	 }
	 else if(n = compArgStr(p+1, "ROTEXOH", 7)){
	    RotExistingOH = TRUE;
	 }
	 else if(n = compArgStr(p+1, "FLIPs", 4)){
	    DemandFlipAllHNQs = TRUE;
	 }
	 else if(n = compArgStr(p+1, "SEGIDmap", 5)){
	    if (++i < argc) {
	       UseSEGIDtoChainMap = TRUE;
	       PDBrec::InstallMapOfSEGIDstoChains(argv[i]);
	    }
	    else {
	       cerr << "no mapping info after -SEGIDmap flag" << endl;
	    }
	 }
	 else if(n = compArgStr(p+1, "Nterm", 1)){
	    MinNTermResNo = parseInteger(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "Model", 1)){
	    ModelToProcess = parseInteger(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "ONLTA", 5)){
	    DoOnlyAltA = TRUE;
	 }
	 else if(n = compArgStr(p+1, "ALLALT", 6)){
	    DoOnlyAltA = FALSE;
	 }
	 else if(n = compArgStr(p+1, "NOHETh", 5)){
	    ProcessConnHydOnHets = FALSE;
	 }
	 else if(n = compArgStr(p+1, "DENSity", 4)){
	    VdwDotDensity = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "PENalty", 3)){
	    PenaltyMagnitude = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "RADius", 3)){
	    ProbeRadius = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "NBonds", 2)){
	    NBondCutoff = parseInteger(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "OCCcutoff", 3)){
	    OccupancyCutoff = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "H2OBcutoff", 4)){
	    WaterBcutoff = 1.0 * parseInteger(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "H2OOCCcutoff", 6)){
	    WaterOCCcutoff = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "HBREGcutoff", 5)){
	    MinRegHBgap = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "HBCHargedcutoff", 4)){
	    MinChargedHBgap = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "BADBumpcutoff", 4)){
	    BadBumpGapCut = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "NONMETALBump", 9)){
	    NonMetalBumpBias = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "METALBump", 6)){
	    MetalBumpBias = parseReal(p, n+1, 10);
	 }
	 else if(n = compArgStr(p+1, "REFerence", 3)){
	    cerr << "Please cite: " << referenceString << endl;
	    cerr << "For more information see " << electronicReference << endl;
	    exit(1);
	 }
	 else if(n = compArgStr(p+1, "FIX", 3)){
	    if (++i < argc) {
			OFile = argv[i];
	    }
	    else {
	       cerr << "no filename after -FIX flag" << endl;
	    }
	 }
	 else if(n = compArgStr(p+1, "DB", 2)){
	    if (++i < argc) {
	       DBfilename = argv[i];
	    }
	    else {
	       cerr << "no filename after -DB flag" << endl;
	    }
	 }
	 else if(n = compArgStr(p+1, "LIMITsearch", 5)){
	    ExhaustiveLimit = parseInteger(p, n+1, 10);
	 }
	 else if(compArgStr(p+1, "Help", 1)){ // has to be after all the other -HXXXs
	    reduceHelp(TRUE);
	 }
	 else {
	    cerr << "unrecognized flag, \"" << p << "\", ignored." << endl;
	 }
      }
      else if (nfile <= 0) {
	    pdbfile = p;
	    nfile = 1;
      }
      else {
	 cerr << "unrecognized parameter, \"" << p << "\", ignored." << endl;
      }
   }
   if (nfile != 1) { reduceHelp(FALSE); }
   return pdbfile;
}

void reduceHelp(bool showAll) { /*help*/
   cerr << versionString << endl;
   cerr << "arguments: [-flags] filename or -" << endl;
   cerr << "040509 reduce.C ln 455: NO renumber, 1459: NO RXR msg" << endl;
   cerr << "041113 rework main to do first and loop over other NMR models if model# not specified."<< endl;
   cerr << "Adds hydrogens to a PDB format file and writes to standard output." << endl;
   cerr << "(note: By default, HIS sidechain NH protons are not added. See -BUILD)" << endl;
   cerr << endl;
   cerr << "Flags:" << endl;
   cerr << "-Trim             remove (rather than add) hydrogens" << endl;
  if (showAll) {
   cerr << endl;
   cerr << "-NOOH             remove hydrogens on OH and SH groups" << endl;
   cerr << "-OH               add hydrogens on OH and SH groups (default)" << endl;
   cerr << endl;
   cerr << "-HIS              create NH hydrogens on HIS rings" << endl;
   cerr << "-FLIPs            allow complete ASN, GLN and HIS sidechains to flip" << endl;
   cerr << "                        (usually used with -HIS)" << endl;
   cerr << "-NOHETh           do not attempt to add NH proton on Het groups" << endl;
   cerr << "-ROTNH3           allow lysine NH3 to rotate (default)" << endl;
   cerr << "-NOROTNH3         do not allow lysine NH3 to rotate" << endl;
   cerr << "-ROTEXist         allow existing rotatable groups (OH, SH, Met-CH3) to rotate" << endl;
   cerr << "-ROTEXOH          allow existing OH & SH groups to rotate" << endl;
   cerr << "-ALLMEthyls       allow all methyl groups to rotate" << endl;
   cerr << "-ONLYA            only adjust 'A' conformations (default)" << endl;
   cerr << "-ALLALT           process adjustments for all conformations" << endl;
   cerr << "-NOROTMET         do not rotate methionine methyl groups" << endl;
   cerr << "-NOADJust         do not process any rot or flip adjustments" << endl;
  }
   cerr << endl;
   cerr << "-BUILD            add H, including His sc NH, then rotate and flip groups" << endl;
   cerr << "                  (except for pre-existing methionine methyl hydrogens)" << endl;
  if (showAll) {
   cerr << "                  (same as: -OH -ROTEXOH -HIS -FLIP)" << endl;
   cerr << endl;
   cerr << "-Keep             keep bond lengths as found" << endl;
   cerr << "-NBonds#          remove dots if cause within n bonds (default="<< NBondCutoff <<")" << endl;
   cerr << "-Model#           which model to process (default="<< ModelToProcess <<")" << endl;
   cerr << "-Nterm#           max number of nterm residue (default="<<MinNTermResNo<<")" << endl;
   cerr << "-DENSity#.#       dot density (in dots/A^2) for VDW calculations (Real, default="<<VdwDotDensity<<")" << endl;
   cerr << "-RADius#.#        probe radius (in A) for VDW calculations (Real, default="<<ProbeRadius<<")" << endl;
   cerr << "-OCCcutoff#.#     occupancy cutoff for adjustments (default="<<OccupancyCutoff<<")" << endl;
   cerr << "-H2OOCCcutoff#.#  occupancy cutoff for water atoms (default="<<WaterOCCcutoff<<")" << endl;
   cerr << "-H2OBcutoff#      B-factor  cutoff for water atoms (Integer, default="<<WaterBcutoff<<")" << endl;
   cerr << "-PENalty#.#       fraction of std. bias towards original orientation (default="<<PenaltyMagnitude<<")" << endl;
   cerr << "-HBREGcutoff#.#   over this gap regular HBonds bump (default="<<MinRegHBgap<<")" << endl;
   cerr << "-HBCHargedcut#.#  over this gap charged HBonds bump (default="<<MinChargedHBgap<<")" << endl;
   cerr << "-BADBumpcut#.#    at this gap a bump is 'bad' (default="<<BadBumpGapCut<<")" << endl;
   cerr << "-METALBump#.#     H 'bumps' metals at radius plus this (default="<<MetalBumpBias<<")" << endl;
   cerr << "-NONMETALBump#.#  'bumps' nonmetal at radius plus this (default="<<NonMetalBumpBias<<")" << endl;
   cerr << "-SEGIDmap \"seg,c...\"  assign chainID based on segment identifier field" << endl;
   cerr << "-Xplor            use Xplor conventions for naming polar hydrogens" << endl;
   cerr << "-NOCon            drop conect records" << endl;
   cerr << "-LIMIT#           max seconds to spend in exhaustive search (default="<< ExhaustiveLimit <<")" << endl;
   cerr << "-NOTICKs          do not display the set orientation ticker during processing" << endl;
   cerr << "-SHOWSCore        display scores for each orientation considered during processing" << endl;
   cerr << "-FIX \"filename\"   if given, file specifies orientations for adjustable groups" << endl;
   cerr << "-DB \"filename\"    file to search for het info" << endl;
   cerr << "                        (default=\""<<DBfilename<<"\")" << endl;
   cerr << "note: can also redirect with unix environment variable: REDUCE_HET_DICT" << endl;
  }
   cerr << endl;
   cerr << "-Quiet            do not write extra info to the console" << endl;
   cerr << "-REFerence        display citation reference" << endl;
   cerr << "-Help             the more extensive description of command line arguments" << endl;
   exit(1);
}

// output a list of PDB records
ostream& outputRecords(ostream& os, const std::list<PDBrec*>& l) {
	for (std::list<PDBrec*>::const_iterator ptr = l.begin(); ptr != l.end(); ++ptr) {
		if ((*ptr)->valid())
			os << (const PDBrec&)(**ptr) << endl;
	}
/*
	DblLnkLstNode<PDBrec>(* ptr)(l._first);  // point to first node in list

  for (int i=0; i < l._num_elem; i++, ptr=l.linkNext(ptr)) {
    if (l.linkData(ptr).valid()) {
      os << (const PDBrec&)(l.linkData(ptr)) << endl;
    }
  }
*/
  return os;
}

// input a list of PDB records
istream& inputRecords(istream& is, std::list<PDBrec*>& records) {
	PDB inputbuffer;
	bool active = TRUE;
	bool modelactive = FALSE;  //041113

	while ((is >> inputbuffer).gcount() > 0) {
		PDBrec* rec = new PDBrec(inputbuffer);
		bool drop = FALSE;

		switch (rec->type()) {
		case PDB::MODEL:
			if (ModelToProcess == rec->modelNum()) {
				active = TRUE;
				modelactive = TRUE; //041113 
			}
			else {
				active = FALSE;
				if(modelactive) //041113
            	{
        	       modelactive = FALSE;
    	           if(ModelSpecified > 0) {ModelNext = 0;} /*only do one*/
	               else {ModelNext = rec->modelNum();} /*next after one just done*/
	            }
				if (Verbose) {
					//cerr << "NOTE: skipping model " << rec->modelNum() << endl;
					if(ModelSpecified > 0)
					{
						cerr << "Model " << ModelToProcess << " specified, skipping model " << rec->modelNum() << endl; //041113
					}
					else
					{
						cerr << "Processing Model " << ModelToProcess << ", for now skipping model " << rec->modelNum() << endl; //041113
					}
			}
		}
		break;
		case PDB::ENDMDL:
			if (! active) { drop = TRUE; active = TRUE; }
			break;
		case PDB::MASTER: drop = TRUE; break; // forget the checksum record
		case PDB::CONECT: if (! KeepConnections) { drop = TRUE; } break;
		case PDB::ATOM: case PDB::HETATM:
			if (active && !drop) {
				if (rec->atomNameModified()) {
					Tally._num_renamed++;
				}
				rec->MapSEGIDtoChain();
			}
			break;
		default:
			break;
		}
		if (active && !drop) {
			records.push_back(rec);
		}
	}
	return is;
}

void renumberAndReconnect(std::list<PDBrec*>& rlst) {
	std::list<PDBrec*>::iterator it = rlst.begin();

	if (KeepConnections && Tally._conect > 0) {
		std::list<PDBrec*> conn;
		std::map<long, PDBrec*> atomsBySeqNum;

		// first we organize atom records by the original sequence num
		// and put all the connect records in a list

		for (; it != rlst.end(); ++it) {
			PDBrec* rin = *it;
			if (rin->type() == PDB::ATOM || rin->type() == PDB::HETATM) {
				if (rin->atomno() > 0) { 
					atomsBySeqNum.insert(std::make_pair(rin->atomno(), rin)); 
				}
			}
			else if (rin->type() == PDB::CONECT) {
				conn.push_back(rin);
			}
		}

		// then we change the numbering

		renumberAtoms(rlst);

		// now we update the connection records

		for (std::list<PDBrec*>::iterator ic = conn.begin(); ic != conn.end(); ++ic) {
			PDBrec* rup = *ic;
			int anum[11];
			rup->getConect(anum);

			for (int k=0; k < 11; k++) {
				if (anum[k] > 0) {
					std::map<long, PDBrec*>::const_iterator iter = atomsBySeqNum.find(anum[k]);
					if (iter != atomsBySeqNum.end())
						anum[k] = iter->second->atomno();
					else 
						anum[k] = -999;
					//	       PDBrec* x = atomsBySeqNum.get(anum[k]);
					//	       anum[k] = (x == NULL) ? -999 : x->atomno();
				}
			}

			rup->setConect(anum);
			if (anum[0] == -999) { rup->invalidateRecord(); }
		}
	}
	else {
		renumberAtoms(rlst);
	}
}

void dropHydrogens(std::list<PDBrec*>& rlst) {
	for (std::list<PDBrec*>::iterator it = rlst.begin(); it != rlst.end(); ++it) {
		PDBrec* r = *it;
		if (r->type() == PDB::ATOM) {
			if (r->isHydrogen()) {
				Tally._H_removed++;
				delete r;
				rlst.erase(it);
				--it;
			}
			Tally._num_atoms++;
		}
		else if (r->type() == PDB::HETATM) {
			if (r->isHydrogen()) {
				Tally._H_removed++;
				Tally._H_HET_removed++;
				delete r;
				rlst.erase(it);
				--it;
			}
			Tally._num_atoms++;
		}
		else if ( r->type() == PDB::SIGATM
			|| r->type() == PDB::ANISOU
            || r->type() == PDB::SIGUIJ) { // supplemental records
			if (r->isHydrogen()) {
				delete r;
				rlst.erase(it);
				--it;
			}
		}
		else if ( r->type() == PDB::CONECT) {
			Tally._conect++;
		}
	}
}

// count record types and group records by positon
// also, update iterator to point to the start of the atom records

void scanAndGroupRecords(std::list<PDBrec*>& rlst, AtomPositions& xyz,
						 std::list<PDBrec*>::iterator& startAtoms) {
	bool foundStart = FALSE;

	for (std::list<PDBrec*>::iterator it = rlst.begin(); it != rlst.end(); ++it) {
		PDBrec* r = *it;
		if (r->type() == PDB::ATOM) {
			if (! foundStart) { foundStart = TRUE; startAtoms = it; }
			if (r->isHydrogen()) {
				Tally._H_found++;
			}
			Tally._num_atoms++;
			// create a new object to avoid delete NULL;
			xyz.put(r);
		}
		else if (r->type() == PDB::HETATM) {
			if (! foundStart) { foundStart = TRUE; startAtoms = it; }
			if (r->isHydrogen()) {
				Tally._H_found++;
				Tally._H_HET_found++;
			}
			Tally._num_atoms++;
			// create a new object to avoid delete NULL;
			xyz.put(r);
		}
		else if (r->type() == PDB::MODEL) {
			if (! foundStart) { foundStart = TRUE; startAtoms = it; }
		}
		else if (r->type() == PDB::CONECT) {
			Tally._conect++;
		}
	}
}

void reduceList(CTab& hetdatabase, std::list<PDBrec*>& rlst,
				AtomPositions& xyz, std::vector<std::string>& fixNotes) {
	std::list<PDBrec*>::iterator it = rlst.begin();
	ResBlk *pb = NULL, *cb = NULL, *nb = NULL;
	std::list<PDBrec*> waters;

	while (it != rlst.end()) { // the ResBlk constructor will advance the iterator
		if (pb) { delete pb; }
		pb = cb;
		cb = nb;

		nb = new ResBlk(rlst, it);
		
		if (nb) { fixupHeavyAtomElementType(*nb, hetdatabase); }

		if (cb && cb->valid(rlst)) {
			analyzeRes(hetdatabase,
				(isConnected(pb, cb) ? pb : NULL),
				cb,
				(isConnected(cb, nb) ? nb : NULL),
				xyz, waters, fixNotes, rlst);
		}
	}
	if (nb && nb->valid(rlst)) {
		analyzeRes(hetdatabase,
			(isConnected(cb, nb) ? cb : NULL), nb, NULL,
			xyz, waters, fixNotes, rlst);
	}

	if (pb) { delete pb; pb = NULL; }
	if (cb) { delete cb; cb = NULL; }
	if (nb) { delete nb; nb = NULL; }

	xyz.generateWaterPhantomHs(waters);
}

void renumberAtoms(std::list<PDBrec*>& rlst) {
	int ia = 0;

	for (std::list<PDBrec*>::iterator it = rlst.begin(); it != rlst.end(); ++it) {
		PDBrec* r = *it;
		if (r->type() == PDB::ATOM
			|| r->type() == PDB::HETATM) {
			r->atomno(++ia);
		}
		else if (r->type() == PDB::SIGATM
            || r->type() == PDB::ANISOU
            || r->type() == PDB::SIGUIJ) {
			r->atomno(ia);
		}
		else if (r->type() == PDB::TER) {
			r->terAtomno(ia);
		}
	}
}

// is the Cterm of a connected to the Nterm of b?
bool isConnected(ResBlk* a, ResBlk* b) {
   if (a == NULL || b == NULL) { return FALSE; }

   std::list<PDBrec*> ar;
   a->get(" C", ar);
   std::list<PDBrec*> br;
   b->get(" N", br);

   if (ar.size() > 0 && br.size() > 0) {

      // must be in same chain
      if ((*(ar.begin()))->chain() != (*(br.begin()))->chain()) { return FALSE; }

      // cterm C of 'a' must be close to nterm N of 'b'
      // we only look at the first set of conformations

      double gap = distance2((*(ar.begin()))->loc(), (*(br.begin()))->loc());
      if (1.1 < gap && gap < 1.7) { return TRUE; } 
   }
   return FALSE;
}

void analyzeRes(CTab& hetdatabase, ResBlk* pb, ResBlk* cb, ResBlk* nb,
				AtomPositions& xyz, std::list<PDBrec*>& waters, 
				std::vector<std::string>& fixNotes, std::list<PDBrec*>& rlst) {
	if (! (cb && cb->valid(rlst))) { return; } // double check

	// add in the previous and next records

	PDBrec* rec_temp = NULL;
	if (pb) {
		std::list<PDBrec*> pr_list;
		pb->get(" C", pr_list);
		for (std::list<PDBrec*>::iterator pr = pr_list.begin(); pr != pr_list.end(); ++pr) {
			cb->addPrevRec(*pr);
		}
	}
	if (nb) {
		std::list<PDBrec*> nr_list;
		nb->get(" N", nr_list);
		for (std::list<PDBrec*>::iterator nr = nr_list.begin(); nr != nr_list.end(); ++nr) {
			cb->addNextRec(*nr);
		}
	}

	const PDBrec rec(cb->firstRec());

	// heuristic to determine if Nterm end of chain or fragment
	ConnType ctype;
	if (pb) { ctype = CONNECTED_RES; }
	else {
		ctype = (rec.resno() <= MinNTermResNo) ? NTERM_RES : FRAGMENT_RES;
	}

	// standardize residue names
	std::string resname = rec.isWater() ? "HOH" : toUppercase(rec.resname());

	if (resname == "HOH") { // add dummy donor protons on waters
		noteWaterInfo(*cb, waters);
	}

	// build flip descriptions
	std::list<char> resAlts;
	if (DemandFlipAllHNQs) { // resAlts will be non-empty for flipped residues
		resAlts = xyz.insertFlip(*cb);
	}

	if (rec.hasProp(METALIC_ATOM)) { // save info from metals for analysis of flips&rots
		xyz.manageMetals(*cb);
	}

	// create plans
	std::list<atomPlacementPlan*> app;

	if (ProcessConnHydOnHets || StdResH::ResXtraInfo().match(resname, "AminoAcid")) {
		StdResH *hn = NULL;
		if (ctype == CONNECTED_RES) {
			hn = StdResH::HydPlanTbl().get("amide");
		}
		else if (ctype == NTERM_RES) {
			bool ispro = resname.find("PRO") != -1;
			hn = StdResH::HydPlanTbl().get(ispro?"nt-pro":"nt-amide");
		}
		else if (ctype == FRAGMENT_RES) { // don't create NH
			if (StandardizeRHBondLengths) {
				// case A - really is a N terminal
				bool ispro = resname.find("PRO") != -1;
				findAndStandardize(ispro?"nt-pro":"nt-amide", resname.c_str(), *cb);

				// case B - really just a fragment
				findAndStandardize("amide", resname.c_str(), *cb);
			}
			// since we don't know for sure if this is an N-terminal or a frag
			// we make sure this nitrogen is an acceptor just in case
			std::list<PDBrec*> nfr_list;
			cb->get(" N", nfr_list);
			for (std::list<PDBrec*>::iterator nfr = nfr_list.begin(); nfr != nfr_list.end(); ++nfr) {
				(*nfr)->elem(* ElementInfo::StdElemTbl().element("Nacc"));
			}
		}
		if (hn && hn->exclude().find(resname.c_str()) == -1) {
			std::list<atomPlacementPlan*> temp = hn->plans();
			app.splice(app.end(), temp);
		}
		if (StdResH::ResXtraInfo().match(resname, "AminoAcid")) {
			StdResH *ha = StdResH::HydPlanTbl().get("alpha");
			if (ha && ha->exclude().find(resname.c_str()) == -1) {
				std::list<atomPlacementPlan*> temp = ha->plans();
				app.splice(app.end(), temp);
			}
		}
	}

	if (StdResH::ResXtraInfo().match(resname, "NucleicAcid")) {
		StdResH *rpbb = StdResH::HydPlanTbl().get("ribose phosphate backbone");
		if (rpbb && rpbb->exclude().find(resname.c_str()) == -1) {
			std::list<atomPlacementPlan*> temp = rpbb->plans();
			app.splice(app.end(), temp);
		}
	}
	bool o2prime = cb->contains(" O2*") || cb->contains(" O2'"); // distinguish RNA and DNA

	// first look in the standard H table...

	//std::cerr << "resname: " << resname << std::endl;
	StdResH *srh = StdResH::HydPlanTbl().get(resname);
	if (srh) {
		//std::cerr << "Found srh for " << resname << std::endl;
		std::list<atomPlacementPlan*> temp = srh->plans();
		app.splice(app.end(), temp);
	}

	else if (ProcessConnHydOnHets) { // if not in the std table we look in the het database...

		ResConn *ct = hetdatabase.findTable(resname);
		if (ct) {
			std::list<atomPlacementPlan*> temp = ct->genHplans();
			app.splice(app.end(), temp);
		}
		else {
			cerr << "*WARNING*: Res \""
				<< resname << "\" not in HETATM Connection Database. Hydrogens not added."
				<< endl;
		}
	}

	// work through each atom placement plan
	for(std::list<atomPlacementPlan*>::const_iterator iter = app.begin(); iter != app.end(); ++iter) {
		genHydrogens(**iter, *cb, o2prime, xyz, resAlts, fixNotes, rlst);
	}
}

void genHydrogens(const atomPlacementPlan& pp, ResBlk& theRes, bool o2prime,
				  AtomPositions& xyz, std::list<char>& resAlts,
				  std::vector<std::string>& fixNotes, std::list<PDBrec*>& rlst) {

	std::list<PDBrec*> ourHydrogens;
	theRes.get(pp.name(), ourHydrogens);
	std::list<PDBrec*> firstAtoms;
	theRes.get(pp.conn(0), firstAtoms);

	bool doNotAdjustSC = FALSE;

	if (!firstAtoms.empty()) { // must connect to something!

		if (!ourHydrogens.empty()) { // Hydrogen exists

			PDBrec* o = NULL;
			for (std::list<PDBrec*>::iterator it = ourHydrogens.begin(); it != ourHydrogens.end(); ++it) {
				o = *it;
				doNotAdjustSC = FALSE;

				if (pp.hasFeature(NOO2PRIMEFLAG) && o2prime) {
					// skip pp record (should find another with the same name)
				}
				else if (pp.hasFeature(O2PRIMEFLAG) && !o2prime) {
					// skip pp record (should find another with the same name)
				}
				else if (pp.hasFeature(UNSUREDROPFLAG)
	                   && ! SaveOHetcHydrogens) {
					noteRemovedInTally(*o);
					o->invalidateRecord();
				}
				else if ( (pp.hasFeature(ROTATEFLAG) &&
					(DemandRotExisting || RotExistingOH))
					|| (pp.hasFeature(ROTATEONDEMAND) &&
					(DemandRotExisting || DemandRotNH3)) ) {

					char hac = o->alt();
					PDBrec r0atom, r1atom, r2atom;        // find connecting atoms
					int connatomcount = findConnAtoms(pp, theRes, hac,
						r0atom, r1atom, r2atom);

					if (connatomcount == 3) {
						xyz.insertRot(*o, r0atom, r1atom, r2atom,
							(DemandRotExisting || RotExistingOH),
							DemandRotNH3,
							(  (DemandRotExisting && DemandRotAllMethyls
							&& pp.hasFeature(ROTATEONDEMAND))
							|| (DemandRotExisting && OKProcessMetMe
							&& pp.hasFeature(ROTATEFLAG))) );
					}

					if (StandardizeRHBondLengths) {
						stdBondLen(pp.dist(), *o,
							firstAtoms, pp.elem());
					}
					if ((connatomcount > 2) &&
						! okToPlaceHydHere(*o, pp, r0atom, r2atom,
						xyz, doNotAdjustSC, fixNotes)) {
						xyz.doNotAdjust(*o);
						noteRemovedInTally(*o);
						o->invalidateRecord();
					}
					else if (doNotAdjustSC) {
						xyz.doNotAdjust(*o);
					}
				}
				else { // non-rotatable hydrogens (incl. flips)
					char hac = o->alt();
					PDBrec r0atom, r1atom, r2atom;        // find connecting atoms
					int connatomcount = findConnAtoms(pp, theRes, hac,
						r0atom, r1atom, r2atom);

					if (StandardizeRHBondLengths) {
						stdBondLen(pp.dist(), *o,
							firstAtoms, pp.elem());
					}
					if ((connatomcount > 2) &&
						! okToPlaceHydHere(*o, pp, r0atom, r2atom,
						xyz, doNotAdjustSC, fixNotes)) {
						xyz.doNotAdjust(*o);
						noteRemovedInTally(*o);
						o->invalidateRecord();
					}
					else if (doNotAdjustSC) {
						xyz.doNotAdjust(*o);
					}
				}
			}
		}
		else {
			int i = 0, j = 0;

			if (pp.hasFeature(NOO2PRIMEFLAG) && o2prime) {
				return; // if o2* atom then we have RNA not DNA
			}
			if (pp.hasFeature(XTRAFLAG) && ! BuildHisHydrogens) {
				return; // do not add NH to His ring without being asked
			}
			if (pp.hasFeature(UNSUREDROPFLAG) && ! SaveOHetcHydrogens) {
				return; // do not do OH, etc. where we can't place reliably
			}
			if ( (pp.hasFeature(NOTXPLORNAME) &&   UseXplorNames)
				||   (pp.hasFeature(   XPLORNAME) && ! UseXplorNames) ) {
				return; // keep our naming convntions straight
			}

			int numConnAtoms = pp.num_conn();
			int maxalt = 0;
			std::vector<int> nconf;
			nconf.reserve(numConnAtoms);

			std::vector< std::vector<PDBrec*> > rvv;
			rvv.reserve(numConnAtoms);
			bool success = TRUE;

			for (i = 0; i < numConnAtoms; i++) { // get the records and count how many
				std::list<PDBrec*> rs;
				theRes.get(pp.conn(i), rs);
				if (!rs.empty()) {
					nconf[i] = rs.size();
					maxalt = max(maxalt, nconf[i]);
					std::vector<PDBrec*> rvv_v;
					rvv_v.reserve(nconf[i]);
					std::list<PDBrec*>::iterator it_rs = rs.begin();
					for(j=0; j < nconf[i]; ++j, ++it_rs) {
						rvv_v.push_back(*it_rs);
						for(int k=j; k > 0; k--) { // sort by altIds
							if ( toupper(rvv_v[j  ]->alt())
								< toupper(rvv_v[j-1]->alt()) ) {
								swap2(rvv_v[j], rvv_v[j-1]);
							}
						}
					}
					rvv.push_back(rvv_v);
				}
				else { success = FALSE; }
			}

			bool considerNonAlt = FALSE;

			if (pp.hasFeature(STRICTALTFLAG) && (numConnAtoms > 3)
				&& (nconf[0] == 1) && (nconf[1] == 1) && (nconf[2] == 1)) {
				maxalt = 1; // this is an H(alpha) so ignore the fourth (CBeta) alternate conf
				considerNonAlt = TRUE;
			}

			// LIMITATION:
			// the logic to determine alt conf codes does not handle the case were the chain
			// of atoms switches codes
			
			std::vector<Point3d> loc(numConnAtoms);
			for(j=0; success && j < maxalt; j++) { // for each alternate conformation...

				char altId = ' ';
				float occ = (*(firstAtoms.begin()))->occupancy();

				for(i=0; i < numConnAtoms; i++) {
					const PDBrec* cnr = rvv[i][min(j, nconf[i]-1)];
					loc[i] = cnr->loc(); //apl 7/3/06 -- FIXING PUSH_BACK BUG

					if (j <= nconf[i]-1) {
						char abc = cnr->alt();   // get the alt loc char to use
						if (abc != ' ' && altId == ' ') { // (use only the first...)
							altId = abc;
							occ = cnr->occupancy();
						}
					}
				}
				if (considerNonAlt) { altId = ' '; occ = (*(firstAtoms.begin()))->occupancy(); }

				Point3d newHpos;

				if (success && (success = pp.placeH(loc, newHpos))) {

					PDBrec* newHatom = new PDBrec();
					(*(firstAtoms.begin()))->clone(newHatom); // dup and modify connected heavy atom

					newHatom->atomname(pp.name().substr(0,4).c_str()); // restrict name size 

					newHatom->x(newHpos.x());
					newHatom->y(newHpos.y());
					newHatom->z(newHpos.z());

					newHatom->annotation("new");

					newHatom->elem(pp.elem());

					newHatom->atomno(0);
					newHatom->alt(altId);

					newHatom->occupancy(occ);

					newHatom->elemLabel(" H");
					newHatom->chargeLabel("  "); // blank since we don't know

					int chrgflgs = basicChargeState(newHatom->atomname(),
						newHatom->resname(),
						PositiveChargeFlag,
						NegativeChargeFlag,
						NULL);
					if (chrgflgs != 0) {
						if ((chrgflgs & PositiveChargeFlag) != 0) {
							newHatom->setPositive();
							newHatom->chargeLabel(" +");
						}
						if ((chrgflgs & NegativeChargeFlag) != 0) {
							newHatom->setNegative();
							newHatom->chargeLabel(" -");
						}
					}
					
					if (visableAltConf(*newHatom, DoOnlyAltA)){ // add hyd. only if visible

						// look for cyclization or other modifications to rot. groups

						if ((numConnAtoms > 2) &&
							! okToPlaceHydHere(*newHatom, pp,
							*(rvv[0][min(j, nconf[0]-1)]),
							*(rvv[2][min(j, nconf[2]-1)]),
							xyz, doNotAdjustSC, fixNotes)) {
							return;  // don't add hyd.
						}

						theRes.insertNewRec(rlst, newHatom);

						noteAddedInTally(*newHatom);

						xyz.put(newHatom); // index in the xyz table

						if (doNotAdjustSC) { return; } // do not add to the adjustable info

						if ( pp.hasFeature(ROTATEFLAG)
							||  (pp.hasFeature(ROTATEONDEMAND)
							&& (DemandRotAllMethyls || DemandRotNH3) )     ) {
							xyz.insertRot(*newHatom,
								*(rvv[0][min(j, nconf[0]-1)]),
								*(rvv[1][min(j, nconf[1]-1)]),
								*(rvv[2][min(j, nconf[2]-1)]),
								TRUE, DemandRotNH3,
								((DemandRotAllMethyls && pp.hasFeature(ROTATEONDEMAND))
								|| (OKProcessMetMe && pp.hasFeature(ROTATEFLAG))) );
						}
						if (DemandFlipAllHNQs && (!resAlts.empty())) {
							xyz.insertFlip(newHatom, resAlts);
						}
					}
				}
			}
		}
	}
}

void noteAddedInTally(const PDBrec& theH) {
   Tally._H_added++;
   if (theH.type() == PDB::HETATM) {
      Tally._H_HET_added++;
   }
}
void noteStdInTally(const PDBrec& theH) {
   Tally._H_standardized++;
   if (theH.type() == PDB::HETATM) {
      Tally._H_HET_standardized++;
   }
}
void noteRemovedInTally(const PDBrec& theH) {
   Tally._H_removed++;
   if (theH.type() == PDB::HETATM) {
      Tally._H_HET_removed++;
   }
}

int findConnAtoms(const atomPlacementPlan& pp, ResBlk& theRes, char hac,
				  PDBrec& r0atom, PDBrec& r1atom, PDBrec& r2atom) {

	int connatomcount = 0;
	PDBrec* rec = NULL;
	
	if (pp.num_conn() > 0) {
		std::list<PDBrec*> r0_list;
		theRes.get(pp.conn(0), r0_list);
		for (std::list<PDBrec*>::iterator r0 = r0_list.begin(); r0 != r0_list.end(); ++r0) {
			rec = *r0;
			char r0ac = rec->alt();
			if (hac == r0ac || hac == ' ' || r0ac == ' ') {
				r0atom = *rec;
				connatomcount++;
				break;
			}
		}
	}

	if (pp.num_conn() > 1) {
		std::list<PDBrec*> r1_list;
		theRes.get(pp.conn(1), r1_list);
		for (std::list<PDBrec*>::iterator r1 = r1_list.begin(); r1 != r1_list.end(); ++r1) {
			rec = *r1;
			char r1ac = rec->alt();
			if (hac == r1ac || hac == ' ' || r1ac == ' ') {
				r1atom = *rec;
				connatomcount++;
				break;
			}
		}
	}

	if (pp.num_conn() > 2) {
		std::list<PDBrec*> r2_list;
		theRes.get(pp.conn(2), r2_list);
		for (std::list<PDBrec*>::iterator r2 = r2_list.begin(); r2 != r2_list.end(); ++r2) {
			rec = *r2;
			char r2ac = rec->alt();
			if (hac == r2ac || hac == ' ' || r2ac == ' ') {
				r2atom = *rec;
				connatomcount++;
				break;
			}
		}
	}
	
	return connatomcount;
}

bool okToPlaceHydHere(const PDBrec& theHatom, const atomPlacementPlan& pp,
					  const PDBrec& a1, const PDBrec& a2, AtomPositions& xyz,
					  bool& doNotAdjustSC, std::vector<std::string>& fixNotes) {

	std::list<PDBrec*> emptySet;

	const Point3d& theHpos = theHatom.loc();

	if (pp.hasFeature(ROTATEFLAG) || pp.hasFeature(NH3FLAG)) {
		const PDBrec& heavyAtom = a1;
		if (heavyAtom.elem().atno() != 6) { // OH, SH and NH3, but not CH3
			const double halfbondlen = heavyAtom.elem().covRad();
			const double  maxbondlen = halfbondlen
				+ ElementInfo::StdElemTbl().maxCovalentRadius();
			std::list<PDBrec*> nearr_list = xyz.neighbors( heavyAtom.loc(),
				halfbondlen + 0.1,
				maxbondlen  + 0.25);
			std::list<PDBrec*> c2batoms;
			int countOfBonds = 0;
			PDBrec* rec = NULL;
			for (std::list<PDBrec*>::const_iterator nearr = nearr_list.begin(); nearr != nearr_list.end(); ++nearr) {
				rec = *nearr;
				if (interactingConfs(heavyAtom, *rec, DoOnlyAltA)
					&& (! rec->elem().isHydrogen())
					&& (! rec->isWater()) ) {
					const double actual = distance2(heavyAtom.loc(),
						rec->loc());
					const double expected = heavyAtom.elem().covRad()
						+ rec->elem().covRad();
					if ((actual >= (expected - 0.55))
						&&  (actual <= (expected + 0.25))) {

						c2batoms.push_front(rec);

						if (++countOfBonds >= 2) {
							// we have an R-X-R or R-X-X-R bond

							if (pp.hasFeature(NH3FLAG)) {
								bool skipthisH = FALSE;
								PDBrec* nbs = NULL;
								for(std::list<PDBrec*>::const_iterator it = c2batoms.begin(); it != c2batoms.end(); ++it) {
									nbs = *it;
									const double hRdist = distance2(theHpos,
										nbs->loc());
									const double bumpLim = nbs->elem().explRad()
										+ NonMetalBumpBias;
									if (hRdist < bumpLim) { skipthisH = TRUE; }

									if ((bumpLim - hRdist) < 0.02) {
										cerr << "WARNING:" << heavyAtom.recName()
											<< " H atom too close to "
											<< nbs->recName() << " by "
											<< (bumpLim - hRdist) << "A" << endl;
										cerr << "        if warning inappropriate, adjust the cuttoffs using "
											<< "-METALBump or -NONMETALBump." << endl;
									}
								}
								
								if (skipthisH) {
									recordSkipInfo(TRUE, fixNotes, theHatom, heavyAtom, c2batoms, "(NH2R)");
									return FALSE; // don't add hyd.
								}
							}
							else { // C-S-S-C, C-O-C, C-O-P, etc.
								//recordSkipInfo(TRUE, fixNotes, theHatom, heavyAtom, c2batoms, "(RXR or RXXR)"); //2.21.mod_dcr
								return FALSE; // don't add hyd.
							}
						}
					}
				}
			}
		}
	}

	// *** special case for N/Q metal coordination ***
	if (pp.hasFeature(BONDBUMPFLAG)
		&& strstr(":ASN:asn:GLN:gln:ASX:asx:GLX:glx:", theHatom.resname())) {
		const PDBrec& nqoxygen = a2;
		const double halfbondlen = nqoxygen.elem().covRad();
		const double  maxbondlen = halfbondlen
			+ ElementInfo::StdElemTbl().maxCovalentRadius();

		std::list<PDBrec*> nearr_list = xyz.neighbors( nqoxygen.loc(),
			halfbondlen + 0.1, maxbondlen  + 0.25);

		PDBrec* rec = NULL;
		for (std::list<PDBrec*>::const_iterator nearr = nearr_list.begin(); nearr != nearr_list.end(); ++nearr) {
			rec = *nearr;
			if (interactingConfs(nqoxygen, *rec, DoOnlyAltA)
				&& rec->hasProp(METALIC_ATOM) ) {
				const double actual = distance2(nqoxygen.loc(), rec->loc());
				const double expected = nqoxygen.elem().covRad() + rec->elem().covRad();

				if ((actual >= (expected - 0.65))
					&& (actual <= (expected + 0.25)) ) {
					// we have determined that the Oxygen is "bonded" to a metal

					emptySet.push_front(rec);
					recordSkipInfo(FALSE, fixNotes, theHatom, nqoxygen, emptySet, "(metal ligand)");
					emptySet.pop_front();
					doNotAdjustSC = TRUE;
				}
			}
		}
	} // *** end of special case code for N/Q metal coord. ***

	if (pp.hasFeature(BONDBUMPFLAG) && (!doNotAdjustSC)) {
		const PDBrec& heavyAtom = a1;
		const double halfbondlen = heavyAtom.elem().covRad();
		const double  maxbondlen = halfbondlen
			+ ElementInfo::StdElemTbl().maxCovalentRadius();

		std::list<PDBrec*> nearr_list = xyz.neighbors( heavyAtom.loc(),
			halfbondlen + 0.1,
			maxbondlen  + 0.25);
		PDBrec* rec = NULL;
		for (std::list<PDBrec*>::const_iterator nearr = nearr_list.begin(); nearr != nearr_list.end(); ++nearr) {
			rec = *nearr;
			if (interactingConfs(heavyAtom, *rec, DoOnlyAltA)
				&& (! (sameres(theHatom, *rec)
				&& StdResH::ResXtraInfo().match(theHatom.resname(), "AminoAcid") ))
				&& (! rec->elem().isHydrogen())
				&& (! rec->isWater()) ) {
				const double actual = distance2(heavyAtom.loc(), rec->loc());
				const double expected = heavyAtom.elem().covRad() + rec->elem().covRad();

				if (interactingConfs(theHatom, *rec, DoOnlyAltA)
					&& (actual >= (expected - 0.65))
					&& (actual <= (expected + 0.25)) ) {

					// this section determines how "bonded" groups place their hydrogens
					// especially as regards how metals determine the protonation pattern

					const double bias = rec->hasProp(METALIC_ATOM)
						? MetalBumpBias     // compensate for small "Metal VDW"
						: NonMetalBumpBias; // skip H very near non-Metal

					const double hRdist = distance2(theHpos, rec->loc());
					const double bumpLim = rec->elem().explRad() + bias;

					if (hRdist < bumpLim) {
						const char *msg = "(H bumps)";
						if ((bumpLim - hRdist) < 0.02) {
							cerr << "WARNING:" << heavyAtom.recName()
								<< " H atom too close to "
								<< rec->recName() << " by "
								<< (bumpLim - hRdist) << "A" << endl;
							cerr << "        you may need to adjust using "
								<< "-METALBump or -NONMETALBump." << endl;
						}
						if (pp.hasFeature(ROTATEFLAG)) { msg = "(short bond)"; }
						if ( rec->hasProp(METALIC_ATOM)
							&& strstr(":ASN:asn:GLN:gln:ASX:asx:GLX:glx:", theHatom.resname())
							&& (heavyAtom.elem().atno()==7)) { // N/Q Nitrogen metal ligand?
							cerr << "WARNING:" << heavyAtom.recName()
								<< " coordinates " << rec->recName()
								<< " (should consider flipping)" << endl;
						}
						emptySet.push_front(rec);
						recordSkipInfo(TRUE, fixNotes, theHatom, heavyAtom, emptySet, msg);
						emptySet.pop_front();
						return FALSE; // don't add hyd.
					}
				}
			}
		}
	}

	if (pp.hasFeature(IFNOPO4)) {
		const PDBrec& theoxygen = a1;
		double pobondlen = theoxygen.elem().covRad() +
			ElementInfo::StdElemTbl().element("P")->covRad();
		std::list<PDBrec*> nearr_list = xyz.neighbors( theoxygen.loc(),
			pobondlen - 0.25, pobondlen + 0.25);
		PDBrec* rec = NULL;
		for (std::list<PDBrec*>::const_iterator nearr = nearr_list.begin(); nearr != nearr_list.end(); ++nearr) {
			rec = *nearr;
			if (interactingConfs(theoxygen, *rec, DoOnlyAltA)
				&& rec->elem().atno() == 15 ) {
				return FALSE; // we have a OP bond -- don't add hyd.
			}
		}
	}
	return TRUE;
}

void recordSkipInfo(bool skipH, std::vector<std::string>& fixNotes,
					const PDBrec& theHatom, const PDBrec& heavyAtom,
					std::list<PDBrec*>& nearr, const char * msg) {
	std::list<PDBrec*> conAtms;
	PDBrec* rec = NULL;
	for(std::list<PDBrec*>::const_iterator it = nearr.begin(); it != nearr.end(); ++it) {
		rec = *it;
		if (! (sameres(heavyAtom, *rec)
			&& StdResH::ResXtraInfo().match(heavyAtom.resname(), "AminoAcid")) ) {
			conAtms.push_front(rec);
		}
	}
	std::string theOtherGuy = (!conAtms.empty()) ? ( (**(conAtms.begin())).recName()
		+ ( (conAtms.size() == 1) ? ":" : ":...") )
		: " cyclic :" ;
	
	std::string heavy = heavyAtom.atomname();

	if ((heavy == " O3'") || (heavy == " O3*") ||
		(heavy == " O5'") || (heavy == " O5*")) {
		for(std::list<PDBrec*>::const_iterator cat = conAtms.begin(); cat != conAtms.end(); ++cat) {
			std::string cat_elem = (*cat)->elemLabel();
			if (cat_elem == " P") {
				return; /* bypass -- do not warn about normal nucleic acid linkages */
			}
		}
	}

	fixNotes.push_back( (skipH ? "USER  MOD NoAdj-H:" : "USER  MOD NoAdj  :") +
		theHatom.recName() + ":" +
		heavyAtom.recName() + ":" +
		theOtherGuy + msg);

	if (Verbose) {
		if (skipH) {
			cerr << "SKIPPED H("<< theHatom.recName() <<"):"
				<< heavyAtom.recName() << " bonds";
		}
		else {
			cerr << "FIXED H("<< theHatom.recName() <<"):"
				<< heavyAtom.recName() << " bonds";
		}
		for(std::list<PDBrec*>::const_iterator it = conAtms.begin(); it != conAtms.end(); ++it) {
			cerr << "-" << (*it)->recName();
		}
		cerr << " " << msg << endl;
	}
}

// For water atoms: ignore any atoms with high b or low occupancy
// and identify waters that do not have explicit H atoms.
// Possible orientations for the new H atoms are identified elsewhere.

void noteWaterInfo(ResBlk& r, std::list<PDBrec*>& waters) {
	std::multimap<std::string, PDBrec*> it_map = r.atomIt();
	std::multimap<std::string, PDBrec*>::const_iterator it = it_map.begin();
	std::string key;
	int i=0, nac=0;
	float occ = 0.0, bf = 0.0;
	const int altbufsz = 10;
	char skipalt[altbufsz], ac;
	bool hwithsamecode;
	PDBrec* a;

	// first gather altcodes for high occupancy Hydrogens
	while(it != it_map.end()) {
		key = it->first;
		if (SaveOHetcHydrogens
			&& (key[0] == ' ' || key[0] == '1' || key[0] == '2')
			&& (key[1] == 'H' || key[1] == 'D')) { // i.e., "[ 12][HD]*"
			for (; it != it_map.end() && it->first == key; ++it) {
				a = it->second;
				occ = a->occupancy();
				ac  = a->alt();
				if (WaterOCCcutoff <= abs(occ)) {
					hwithsamecode = FALSE;
					for(i=0; i < nac; i++) {
						if (skipalt[i] == ac) { hwithsamecode = TRUE; break; }
					}
					if ((! hwithsamecode) && (nac < altbufsz)) {
						skipalt[nac++] = ac;
					}
				}
			}
		}
		else
			++it;
	}

	// then, in a second pass, we cull out the low occ and high b atoms
	// and for the remaining oxygens note if there are no
	// high occupancy protons with the same alt code.
	it = it_map.begin();
	while(it != it_map.end()) {
		key = it->first;
		for (; it != it_map.end() && it->first == key; ++it) {
			a = it->second;
			occ = a->occupancy();
			bf  = a->tempFactor();
			if (WaterOCCcutoff <= abs(occ) && bf < WaterBcutoff) {
				if (key[0] == ' ' && key[1] == 'O') { // i.e., " O*"
					hwithsamecode = FALSE;
					for(i=0; i < nac; i++) {
						if (skipalt[i] == a->alt()) { // requires exact match
							hwithsamecode = TRUE;
							break;
						}
					}
					if (! hwithsamecode) { waters.push_front(a); } //save for later...
				}
			}
			else { // ignore water atoms which are not good enough
				a->elem(* ElementInfo::StdElemTbl().element("ignore"));
			}
		} // end for loop
	}
}

// find Hydrogen Plans by name and standardize bond lengths for those Hs
void findAndStandardize(const char* name, const char* resname,
						ResBlk& theRes) {
	StdResH *hptr = StdResH::HydPlanTbl().get(name);

	if (hptr && hptr->exclude().find(resname) == -1) {
		std::list<atomPlacementPlan*> app_deque = hptr->plans();
		for (std::list<atomPlacementPlan*>::iterator app = app_deque.begin(); app != app_deque.end(); ++app) {

			const atomPlacementPlan* pp = *app;
			std::list<PDBrec*> firstAtoms;
			theRes.get(pp->conn(0), firstAtoms);

			if (!firstAtoms.empty()) { // must connect to something!

				std::list<PDBrec*> ourHydrogens;
				theRes.get(pp->name(), ourHydrogens);
				for (std::list<PDBrec*>::iterator it = ourHydrogens.begin(); it != ourHydrogens.end(); ++it) {
					stdBondLen(pp->dist(), **it, firstAtoms, pp->elem());
				}
			}
		}
	}
}

void stdBondLen(float dist, PDBrec& ourHydrogen, std::list<PDBrec*>& firstAtoms,
				const ElementInfo& e) {
	if (ourHydrogen.valid()) {
		PDBrec* temp = NULL;
		for(std::list<PDBrec*>::iterator it = firstAtoms.begin(); it != firstAtoms.end(); ++it) {
			temp = *it;
			char halt = ourHydrogen.alt();
			char calt = temp->alt();
			if (halt == calt || halt == ' ' || calt == ' ') {
				Point3d hpos = ourHydrogen.loc();
				Point3d cpos = temp->loc();

				Vector3d bondvec = hpos - cpos;
				Point3d  newHpos = cpos + bondvec.scaleTo(dist);

				ourHydrogen.x(newHpos.x());
				ourHydrogen.y(newHpos.y());
				ourHydrogen.z(newHpos.z());
				
				float bfactor = ourHydrogen.tempFactor();
				if (bfactor < 1.0) {     // inherit bfactor if none given
					ourHydrogen.tempFactor(temp->tempFactor());
				}

				ourHydrogen.annotation("std");

				if (! ourHydrogen.hasProp(IGNORE)) {
					ourHydrogen.elem(e);             // specialize element type
				}

				noteStdInTally(ourHydrogen);
			}
		}
	}
}

void fixupHeavyAtomElementType(ResBlk& theRes, CTab& hetDB) {
	std::multimap<std::string, PDBrec*> theAtoms_map = theRes.atomIt();
	std::multimap<std::string, PDBrec*>::const_iterator theAtoms = theAtoms_map.begin();
	std::string atomname;
	while (theAtoms != theAtoms_map.end()) {
		atomname = theAtoms->first;
		PDBrec* r;
		for (; theAtoms != theAtoms_map.end() && theAtoms->first == atomname; ++theAtoms) {
			r = theAtoms->second;

			// Currently, all we are fixing up is the HBonding status of specific carbon & nitrogen atoms.
			// Elsewhere we add an expanded HOd donor atom around the oxygen of water to provide
			// an omnidirectional Hbonding arrangement
			if (r->hasProp(IGNORE)) {}          // don't alter these
			else if (r->elem().atno() == 7) {                              // nitrogen
				if (StdResH::ResXtraInfo().atomHasAttrib(r->resname(), atomname, HBACCEPTORFLAG)) {
					r->elem(* ElementInfo::StdElemTbl().element("Nacc"));
				}
				else {
					const int nc = hetDB.numConn(r->resname(), atomname);
					if (nc == 1 || nc == 2) {            // HET NR or NR2 H-Bond acceptor
						r->elem(* ElementInfo::StdElemTbl().element("Nacc"));
					}
				}
			}
			else if (r->elem().atno() == 6) {                              // carbon
				if (StdResH::ResXtraInfo().atomHasAttrib(r->resname(), atomname, AROMATICFLAG)) {
					r->elem(* ElementInfo::StdElemTbl().element("Car"));
				}
				else if ((atomname == " C")
					&& StdResH::ResXtraInfo().match(r->resname(), "AminoAcid")) {
					r->elem(* ElementInfo::StdElemTbl().element("C=O"));
				}
				else if (StdResH::ResXtraInfo().atomHasAttrib(r->resname(), atomname, ISACOFLAG)) {
					r->elem(* ElementInfo::StdElemTbl().element("C=O"));
				}
				// ***** in addition, aromatic and carbonyl HET carbons need to be identified somehow
			}
		}
	}
}
