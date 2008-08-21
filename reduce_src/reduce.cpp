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
// Copyright (C) 1999-2008 J. Michael Word
// **************************************************************
//
//  reduceChanges now contains the CHANGELOG or history info
//
 
#if defined(_MSC_VER)
#pragma warning(disable:4786) 
#pragma warning(disable:4305) 
#pragma warning(disable:4800) 
#endif

static const char *versionString =
     "reduce: version 3.14 08/21/2008, Copyright 1997-2008, J. Michael Word";

static const char *shortVersion    = "reduce.3.14.080821";
static const char *referenceString =
                       "Word, et. al. (1999) J. Mol. Biol. 285, 1735-1747.";
static const char *electronicReference = "http://kinemage.biochem.duke.edu";

#include <iostream>
#include <fstream>
#include <sstream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#ifdef OLD_STD_HDRS
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#else
#include <cstdlib>
#include <cstring>
#include <cctype>
using std::exit;
using std::strstr;
using std::toupper;
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
bool UseOldNames	      = FALSE; 
bool BackBoneModel	      = FALSE; 
bool DemandRotAllMethyls      = FALSE;
bool RotExistingOH            = FALSE;
bool NeutralTermini           = FALSE;
bool DemandRotNH3             = TRUE;
bool DemandRotExisting        = FALSE;
bool DemandFlipAllHNQs        = FALSE;
bool DoOnlyAltA               = TRUE;
bool OKProcessMetMe           = TRUE;
bool OKtoAdjust               = TRUE;
bool ShowCliqueTicks          = TRUE;
bool ShowOrientScore          = FALSE;
bool StringInput              = FALSE;
bool ShowCharges              = FALSE;

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
float GapWidth        = 0.3; // half width for detecting chain breaks between residues 
                             // (center at 1.4; default allow 1.1-1.7 for accepting connected residues) 

std::string OFile; // if file exists, given orientations forced
bool UseSEGIDtoChainMap = FALSE; // if true, override some chain ids

#ifndef HET_DICTIONARY
//#define HET_DICTIONARY "reduce_het_dict.txt"
#define HET_DICTIONARY "reduce_wwPDB_het_dict.txt"
#endif
#ifndef HET_DICTOLD
#define HET_DICTOLD "reduce_het_dict.txt"
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
void processPDBfile(std::istream& ifs, char *pdbFile, std::ostream& ofs);
void establishHetDictionaryFileName(void);
void reduceHelp(bool showAll);
void reduceChanges(bool showAll);
std::istream& inputRecords(std::istream& is, std::list<PDBrec*>& records);
std::ostream& outputRecords(std::ostream& os, const std::list<PDBrec*>& records);
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
//   std::ifstream theinputstream(pdbFile);
//   processPDBfile(theinputstream, pdbFile, std::ofstream("dump.out"));
//
//   return ReturnCodeGlobal; // one pass and then we quit
//}
//#else

int main(int argc, char **argv) {

   char *pdbFile = parseCommandLine(argc, argv);

   if(pdbFile && StringInput == FALSE)
   {/*can do multiple passes by rewinding input file*/
      //std::ifstream theinputstream(pdbFile); //would need rewind
      while(ModelToProcess) /* 041113 */
      {
         std::ifstream theinputstream(pdbFile); //declare each time, avoid rewind
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
   /* pdbFile here is used to hold the string to be memory efficient */
   else if(StringInput == TRUE){
       std::istringstream is(pdbFile); 
       processPDBfile(is,NULL, cout);
   }
   else
   {/*presume stdin for pdb file info*/
      processPDBfile(cin,            pdbFile, cout); 
   }
   return ReturnCodeGlobal;
}
//#endif

void processPDBfile(std::istream& ifs, char *pdbFile, std::ostream& ofs) {
   if (Verbose) {
      cerr << versionString << endl;
      if (pdbFile) {
	 cerr << "Processing file: \"" << pdbFile << "\"" << endl;
      }
      else if (StringInput){
          cerr << "Processing input string" << endl;
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
            if (NeutralTermini) {
               cerr << "Adding \"amide\" Hydrogens to chain breaks." << endl; 
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

      AtomPositions xyz(2000, DoOnlyAltA, UseXplorNames, UseOldNames, BackBoneModel, 
			NBondCutoff, MinRegHBgap, MinChargedHBgap, 
			BadBumpGapCut, dotBucket, ProbeRadius,
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
	 else if(compArgStr(p+1, "STRING", 6)){
         StringInput = TRUE;
     }
     else if(compArgStr(p+1, "BUILD", 5)){
	    BuildHisHydrogens  = TRUE;
	    SaveOHetcHydrogens = TRUE;
	    RotExistingOH      = TRUE;
	    DemandFlipAllHNQs  = TRUE;
         }
         else if(n = compArgStr(p+1, "NOBUILD", 7)){
            PenaltyMagnitude = parseReal(p, n+1, 10);
         // PenaltyMagnitude = 200;      9999 in molprobity
            BuildHisHydrogens  = TRUE;
            SaveOHetcHydrogens = TRUE;
            RotExistingOH      = TRUE;  //  not used in molprobity
            DemandFlipAllHNQs  = TRUE;
         }
	 else if((n = compArgStr(p+1, "Version", 1))){
	    cerr << shortVersion << endl; 
	    exit(1); 
	 }
	 else if((n = compArgStr(p+1, "Changes", 1))) {
     	     reduceChanges(TRUE);
         }
	 else if((n = compArgStr(p+1, "Quiet", 1))){
	    Verbose = FALSE;
	 }
	 else if((n = compArgStr(p+1, "NOTICKs", 6))){
	    ShowCliqueTicks = FALSE;
	 }
	 else if((n = compArgStr(p+1, "SHOWSCore", 6))){
	    ShowOrientScore = TRUE;
	 }
	 else if((n = compArgStr(p+1, "NOCon", 3))){
	    KeepConnections = FALSE;
	 }
	 else if((n = compArgStr(p+1, "NOROTMET", 8))){
	    OKProcessMetMe = FALSE;
	 }
	 else if((n = compArgStr(p+1, "NOADJust", 5))){
	    OKtoAdjust = FALSE;
	 }
	 else if((n = compArgStr(p+1, "HIS", 3))){
	    BuildHisHydrogens = TRUE;
	 }
	 else if((n = compArgStr(p+1, "OH", 2))){
	    SaveOHetcHydrogens = TRUE;
	 }
	 else if((n = compArgStr(p+1, "NOOH", 4))){
	    SaveOHetcHydrogens = FALSE;
	 }
	 else if((n = compArgStr(p+1, "Xplor", 1)) && !UseOldNames){
	    UseXplorNames = TRUE;
	 }
         else if((n = compArgStr(p+1, "Xplor", 1)) && UseOldNames){
            cerr << "Cannot use both -Xplor and -OLDpdb flags" << endl;
            exit(1);         
	 }
	 else if((n = compArgStr(p+1, "OLDpdb", 3)) && ! UseXplorNames){
	    UseOldNames = TRUE; 
	    DBfilename = HET_DICTOLD;
	 }
         else if((n = compArgStr(p+1, "OLDpdb", 3)) && UseXplorNames){
	    cerr << "Cannot use both -Xplor and -OLDpdb flags" << endl; 
	    exit(1);
	 }
	 else if((n = compArgStr(p+1, "BBmodel", 2))){
		BackBoneModel = TRUE; 
	 } 
	 else if((n = compArgStr(p+1, "Trim", 1))){
	    RemoveHydrogens = TRUE;
	 }
	 else if((n = compArgStr(p+1, "Keep", 1))){
	    StandardizeRHBondLengths = FALSE;
	 }
	 else if((n = compArgStr(p+1, "ALLMETHYLS", 10))){
	    DemandRotAllMethyls = TRUE;
	 }
	 else if((n = compArgStr(p+1, "ROTEXist", 5))){
	    DemandRotExisting = TRUE;
	 }
         else if((n = compArgStr(p+1, "ADDNHATGAP", 10))){
            NeutralTermini = TRUE;
         }
	 else if((n = compArgStr(p+1, "ROTNH3", 6))){
	    DemandRotNH3 = TRUE;
	 }
	 else if((n = compArgStr(p+1, "NOROTNH3", 8))){
	    DemandRotNH3 = FALSE;
	 }
	 else if((n = compArgStr(p+1, "ROTEXOH", 7))){
	    RotExistingOH = TRUE;
	 }
	 else if((n = compArgStr(p+1, "FLIPs", 4))){
	    DemandFlipAllHNQs = TRUE;
	 }
	 else if((n = compArgStr(p+1, "SEGIDmap", 5))){
	    if (++i < argc) {
	       UseSEGIDtoChainMap = TRUE;
	       PDBrec::InstallMapOfSEGIDstoChains(argv[i]);
	    }
	    else {
	       cerr << "no mapping info after -SEGIDmap flag" << endl;
	    }
	 }
	 else if((n = compArgStr(p+1, "Nterm", 1))){
	    MinNTermResNo = parseInteger(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "Model", 1))){
	    ModelToProcess = parseInteger(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "ONLTA", 5))){
	    DoOnlyAltA = TRUE;
	 }
	 else if((n = compArgStr(p+1, "ALLALT", 6))){
	    DoOnlyAltA = FALSE;
	 }
     else if((n = compArgStr(p+1, "CHARGEs", 6))){
        ShowCharges = TRUE;
     }
	 else if((n = compArgStr(p+1, "NOHETh", 5))){
	    ProcessConnHydOnHets = FALSE;
	 }
	 else if((n = compArgStr(p+1, "DENSity", 4))){
	    VdwDotDensity = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "PENalty", 3))){
	    PenaltyMagnitude = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "RADius", 3))){
	    ProbeRadius = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "NBonds", 2))){
	    NBondCutoff = parseInteger(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "OCCcutoff", 3))){
	    OccupancyCutoff = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "H2OBcutoff", 4))){
	    WaterBcutoff = 1.0 * parseInteger(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "H2OOCCcutoff", 6))){
	    WaterOCCcutoff = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "HBREGcutoff", 5))){
	    MinRegHBgap = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "HBCHargedcutoff", 4))){
	    MinChargedHBgap = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "BADBumpcutoff", 4))){
	    BadBumpGapCut = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "NONMETALBump", 9))){
	    NonMetalBumpBias = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "METALBump", 6))){
	    MetalBumpBias = parseReal(p, n+1, 10);
	 }
         else if((n = compArgStr(p+1, "GAPERROR", 8))){
            GapWidth = parseReal(p, n+1, 10);
            if (GapWidth > 1.4) {
               cerr << "Max allowed HalfGapWidth is 1.4" << endl; 
               exit(1); 
            }
         }
	 else if((n = compArgStr(p+1, "REFerence", 3))){
	    cerr << "Please cite: " << referenceString << endl;
	    cerr << "For more information see " << electronicReference << endl;
	    exit(1);
	 }
	 else if((n = compArgStr(p+1, "FIX", 3))){
	    if (++i < argc) {
			OFile = argv[i];
	    }
	    else {
	       cerr << "no filename after -FIX flag" << endl;
	    }
	 }
	 else if((n = compArgStr(p+1, "DB", 2))){
	    if (++i < argc) {
	       DBfilename = argv[i];
	    }
	    else {
	       cerr << "no filename after -DB flag" << endl;
	    }
	 }
	 else if((n = compArgStr(p+1, "LIMITsearch", 5))){
	    ExhaustiveLimit = parseInteger(p, n+1, 10);
	 }
	 else if(compArgStr(p+1, "Help", 1)){ // has to be after all the other -HXXXs
	    reduceHelp(TRUE);
	 }
	 else {
	    cerr << "unrecognized flag, \"" << p << "\", ignored." << endl;
	 }
      }
      else if (StringInput == TRUE){
          pdbfile = p;
          nfile = 1;
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
   cerr << shortVersion << endl; 
   cerr << "arguments: [-flags] filename or -" << endl;
//   cerr << "040509 reduce.C ln 455: NO renumber, 1459: NO RXR msg" << endl;
//   cerr << "041113 rework main to do first and loop over other NMR models if model# not specified."<< endl;
//   cerr << "Adds hydrogens to a PDB format file and writes to standard output." << endl;
//   cerr << "(note: By default, HIS sidechain NH protons are not added. See -BUILD)" << endl;
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
//   cerr << "-ADDNHATGAP            add \"amide\" hydrogen on chain breaks" <<endl; 
//   cerr << "-GAPERROR#.#       sets the half width for allowed peptide bond lengths variations around 1.4 Angstroms: default 0.3" << endl; 
   cerr << "-ROTNH3           allow lysine NH3 to rotate (default)" << endl;
   cerr << "-NOROTNH3         do not allow lysine NH3 to rotate" << endl;
   cerr << "-ROTEXist         allow existing rotatable groups (OH, SH, Met-CH3) to rotate" << endl;
   cerr << "-ROTEXOH          allow existing OH & SH groups to rotate" << endl;
//   cerr << "-ALLMETHYLS       allow all methyl groups to rotate" << endl;
   cerr << "-ONLYA            only adjust 'A' conformations (default)" << endl;
   cerr << "-ALLALT           process adjustments for all conformations" << endl;
   cerr << "-CHARGEs          output charge state for appropriate hydrogen records" << endl;
   cerr << "-NOROTMET         do not rotate methionine methyl groups" << endl;
   cerr << "-NOADJust         do not process any rot or flip adjustments" << endl;
  }
   cerr << endl;
   cerr << "-NOBUILD#.#       build with a given penalty often 200 or 999" << endl;
   cerr << "-BUILD            add H, including His sc NH, then rotate and flip groups" << endl;
   cerr << "                  (except for pre-existing methionine methyl hydrogens)" << endl;
   cerr << endl;
  if (showAll) {
   cerr << "                  (same as: -OH -ROTEXOH -HIS -FLIP)" << endl;
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
   cerr << "-OLDpdb 	      use the pre-remediation names for hydrogens" << endl; 
   cerr << "-BBmodel	      expects a backbone only model and will build HA hydrogens on Calpha truncated residues" <<endl; 
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
   cerr << "-STRING           pass reduce a string object from a script, must be quoted" << endl;
   cerr << "usage: from within a script, reduce -STRING \"_name_of_string_variable_\"" << endl << endl;
   cerr << "-Quiet            do not write extra info to the console" << endl;
   cerr << "-REFerence        display citation reference" << endl;
   cerr << "-Version          display the version of reduce" <<endl; 
   cerr << "-Changes          display the change log" <<endl;
   cerr << "-Help             the more extensive description of command line arguments" << endl;
   exit(1);
}

void reduceChanges(bool showAll) { /*changes*/
   cerr  <<  "History: this new version of reduce has been completely re-written" << endl; 
   cerr  <<  "  to include het groups and rotations" << endl ; 
   cerr  << endl; 
   cerr  << "10/ 1/97 - jmw - changed rotation score to dot score, modified geometry," << endl; 
   cerr  << "                 of CH2 and mcNH angles to be compatible with ECEPP," << endl; 
   cerr  << "                 and added option to permit all methyls to rotate." << endl; 
   cerr  << "10/22/97 - jmw - fixed bug in rotation code which did not separate" << endl; 
   cerr  << "                 alternate conformations" << endl; 
   cerr  << "11/12/97 - jmw - support for flips of HIS and ASN/GLN residues," << endl; 
   cerr  << "                 orientable waters and Hbonds to aromatics" << endl; 
   cerr  << "                 (includes some analysis of pairing of flips and rots)" << endl; 
   cerr  << "12/10/97 - jmw - fixed waters to have phantom H atoms," << endl; 
   cerr  << "                 better rotational searching, metal covalent radii" << endl; 
   cerr  << " 2/ 7/98 - jmw - broke out motions and re-wrote search," << endl; 
   cerr  << "                 use only A conf, fix penalty" << endl; 
   cerr  << " 2/18/98 - jmw - fixup ambiguous atom names, no creation of altB Hs," << endl; 
   cerr  << "                 MTO waters, added -BUILD and other flags," << endl; 
   cerr  << "                 better recognise S-S or C-O-C bonds, etc." << endl; 
   cerr  << " 3/ 1/98 - jmw - fix orientations" << endl; 
   cerr  << " 3/25/98 - jmw - updated metal radii, added HET dict env var," << endl; 
   cerr  << "                 expanded fixed orientations, changed -OH default" << endl; 
   cerr  << " 4/ 8/98 - jmw - fixed bug, occupancy of water phantomHs now > 0.0 !" << endl; 
   cerr  << " 4/18/98 - jmw - bad bump check, short water h, new limit on HB, ..." << endl; 
   cerr  << " 4/26/98 - jmw - new hires search strategy" << endl; 
   cerr  << " 4/29/98 - jmw - fixed min rotation angle and penalty bugs" << endl; 
   cerr  << " 5/13/98 - jmw - fixed another min rotation angle bug where the coarse score" << endl; 
   cerr  << "                 is the best that can be obtained. Also updated output to" << endl; 
   cerr  << "                 the header to document K/C/X/F categories." << endl; 
   cerr  << " 5/15/98 - jmw - fine tuned cases where group is fixed by metal or modification" << endl; 
   cerr  << " 7/ 3/98 - jmw - fixed bug: ASN/GLN sc C=O carbon radii was 1.75, now 1.65" << endl; 
   cerr  << " 8/ 7/98 - jmw - stop putting notes in the segID, elem & charge fields" << endl; 
   cerr  << "                 because PDB use of these fields are now being expanded" << endl; 
   cerr  << "                 by naming convention differences with XPLOR and XPLORs" << endl; 
   cerr  << "                 use of SEGID rather than CHAINID" << endl; 
   cerr  << " 8/12/98 - jmw - changed how we figure num of permutations" << endl; 
   cerr  << " 8/14/98 - jmw - now search for het dictionary in current directory" << endl; 
   cerr  << " 9/ 1/98 - jmw - worked on portability by cleaning up g++ warnings" << endl; 
   cerr  << " 9/13/98 - jmw - figured out a work-around to g++ template and" << endl; 
   cerr  << "                 static const class initializer problems" << endl; 
   cerr  << " 1/ 7/99 - jmw - consolidate Linux and Mac changes with SGI src" << endl; 
   cerr  << " 3/16/99 - jmw - extended atom name parsing for wierd HETs with col1 ABCDEFGs" << endl; 
   cerr  << " 4/ 6/99 - jmw - added OW to list of names for oxygen in water" << endl; 
   cerr  << " 7/ 7/99 - jmw - Improved portability (near->nearr, List <T>::, Makefiles)," << endl; 
   cerr  << "                 made penalty 1.0, added reference, changed basic io hooks" << endl; 
   cerr  << " 9/ 2/99 - jmw - Updated Makefile list and Utility.h for DEC alpha support" << endl; 
   cerr  << "10/ 5/99 - jmw - Updated a Makefiles and the main in reduce.C for sgi6.5 compiler" << endl; 
   cerr  << " 8/ 4/00 - jmw - Added -segid flag to support segment identifiers" << endl; 
   cerr  << "10/30/00 - jmw - Modified Seq mergesort const/non const stuff" << endl; 
   cerr  << " 4/19/01 - jmw - Added support for left justified A/T/U/C/G for nucleic acids" << endl; 
   cerr  << "                 which fixed a rare but nasty bug in the water recognition" << endl; 
   cerr  << " 5/ 7/01 - jmw - Stopped output with -quiet while fixing orientation" << endl; 
   cerr  << "                 and finally dealt with string literal conversion messages" << endl; 
   cerr  << " 5/31/01 - jmw - Pass signal of abandoned clique search in rc," << endl; 
   cerr  << "                 added new flag to allow mapping of segids to chains," << endl; 
   cerr  << "                 changed properties of Nterminal fragment nitrogens to acceptor" << endl; 
   cerr  << "10/ 4/01 - jmw - fiddled to get compiled on RH linux 7 (gcc 2.96)" << endl; 
   cerr  << " 5/24/02 - jmw - added control over hbump" << endl; 
   cerr  << " 1/23/03 - jmw -v2.16 - Changed global MinChargedHBgap: 0.4 => 0.8" << endl; 
   cerr  << "                        Changed phosphorus properties to drop acceptor status" << endl; 
   cerr  << " 3/06/03 - jmw -v2.17 - Cleaning declarations found by Leo and Andrew such as" << endl; 
   cerr  << "                        adding compile flag -DOLD_CXX_DEFNS to select const longs" << endl; 
   cerr  << "                        instead of std::_Ios_Fmtflags in AtomPositions.C" << endl; 
   cerr  << " 3/31/03 - jmw -v2.18 - Fixed spelling mistake by renaming cuttoff to cutoff throughout." << endl; 
   cerr  << "                        Edited help for -H2OBcutoff# to suggest an integer value" << endl; 
   cerr  << "                        Made -H2OOCCcutoff#.# parse a real (was integer)" << endl; 
   cerr  << " 4/ 3/03 - jmw -v2.19 - Moved -Help to the end of the command line parsing." << endl; 
   cerr  << " 6/ 3/03 - jmw -v2.20 - updated to isoC++ style std includes: #include <cstring>," << endl; 
   cerr  << "                        changed String class to Stringclass class," << endl; 
   cerr  << "                        nice speedups from Jack Snoeyink" << endl; 
   cerr  << " 6/ 4/03 - jmw -v2.21 - added out for OLD_STD_HDRS for sgi plus other sgi polishing" << endl; 
   cerr  << "11/ 7/03 - jmw -v2.22 - fixed bug with creating hydrogens at the end of a triple bond" << endl; 
   cerr  << "                        and supressed the warnings about connecting nucleic acid bases" << endl; 
   cerr  << endl; 
   cerr  << "040509dcr reduce.C ln 455: NO renumber, 1459: NO RXR msg" << endl; 
   cerr  << ""<< endl; 
   cerr  << "041113dcr reduce.C reconstructed main to loop over any NMR models in the file" << endl; 
   cerr  << "          a specific model can still be specified." << endl; 
   cerr  << endl; 
   cerr  << "changes" << endl; 
   cerr  << endl; 
   cerr  << "050314 dcr incorporated Mike's version of 030703, to wit:" << endl; 
   cerr  << endl; 
   cerr  << " 6/15/06 - apl -v3.0  - incorporated decomposition of scoring function into interactions of" << endl; 
   cerr  << "                        atom singles, atom pairs, atom tripples and all other" << endl; 
   cerr  << "                        higher-order atom interactions.  incorporated" << endl; 
   cerr  << "                        dynamic programming. incorporated changes from v2.21.mod_dcr" << endl; 
   cerr  << "                        disabling hydrogen atom numbering and skipInfo. removing" << endl; 
   cerr  << "                        PDBrecNAMEout as in v2.999. Disabling hydrogen bonds to his" << endl; 
   cerr  << "                        aromatic carbon atoms." << endl; 
   cerr  << " 6/19/06 - apl -v3.01- decomposing the scoring function in terms of which dots should" << endl; 
   cerr  << "                       be scored with which hyperedges.  Incorporating S3 reduction rules" << endl; 
   cerr  << "                       into dynamic programming.  Incorporating code to handle 4-way overlap" << endl; 
   cerr  << "                       (but not five way overlap) though I have not observed any 4-way" << endl; 
   cerr  << "                       overlap using the new decomposition scheme (it would show up in the previous scheme)." << endl; 
   cerr  << " 6/24/06 - apl -v3.02- incorporating the additions to v2.23 that deal with multiple NMR models" << endl; 
   cerr  << "                       in a single file.  Altering output from dp to be less intrusive and a" << endl; 
   cerr  << "                       little more informative." << endl; 
   cerr  << "                       Adding new reduction rule for vertices with exactly 1 state (i.e. no" << endl; 
   cerr  << "                       real options) which speeds up dynamic programming for the second-round" << endl; 
   cerr  << "                       of optimizations that calculate the optimal network states for sub-optimal" << endl; 
   cerr  << "                       flip states.  GLN and ASN will have 1 state in these optimizations." << endl; 
   cerr  << "7/03/06 - apl -      - Fixing USER MOD Set ordering in output PDB.  Fixing 'flip' records in columns" << endl; 
   cerr  << "                       82 to 85 for the atoms on flipped residues.  Fixing atom placement plan bug" << endl; 
   cerr  << "                       in genHydrogens() that manifested itself as aberant behavior when presented" << endl; 
   cerr  << "                       with several (3 or more) alternate conformations." << endl; 
   cerr  << "7/09/06 - apl -      - Adding #include <cassert> for gcc3.3.4 builds" << endl; 
   cerr  << "7/11/06 - apl -      - Fixing 'node with 0 states' bug." << endl; 
   cerr  << "10/19/06 - apl - v3.03- changing HIS carbons to regular carbons and not arromatic carbons" << endl; 
   cerr  << "10/20/06 - apl -      - fixing bug in optimization code that failed to keep optimal network states" << endl; 
   cerr  << "                       when a network was forced to incur a penalty." << endl; 
   cerr  << "3/ 7/07 - apl -        Bug fix: do not add hydrogens to 'N' on non-amino acids" << endl; 
   cerr  << "                       fixes 3H bug on SAC in 1b0b.pdb" << endl; 
   cerr  << "3/ 7/07 - apl -        Bug fix: march ResBlk::_insertPtr backwards at the end of constructor" << endl; 
   cerr  << "                       even when _insertPtr has reached the end of the rlst.  This ensures all" << endl; 
   cerr  << "                       residues are protonated, instead of all but the last one.  Fixes Loren's" << endl; 
   cerr  << "                       bug in trying to protonate nicotine when it was by itself in a .pdb." << endl; 
   cerr  << "3/ 7/07 - apl -        inline distanceSquared in toolclasses/Point3d.h for a 10% speedup." << endl; 
   cerr  << endl; 
   cerr  << "7/ 7/07 - jjh & rmi -  incoporated new hydrogen names in StdResH.cpp for remediated pdb files" << endl; 
   cerr  << "               v3.10   followed the example of the -Xplor flag and added a -OLDpdb flag to allow" << endl; 
   cerr  << "                       output of new (default) or old (pre-remediation) hydrogen names" << endl; 
   cerr  << "7/13/07 - jjh & rmi -  added -BBmodel flag allows addition of hydrogens on Calpha of truncated" << endl; 
   cerr  << "                       amino acids" << endl; 
   cerr  << "7/16/07 - jjh & rmi -  added fixAtomName() to ElementInfo.cpp to check for Hg, Ho, and Hf atoms" << endl; 
   cerr  << "                       the corresponding warnings are no longer output" << endl; 
   cerr  << "7/31/07 - rmi -        Bug fix: changed the logic in selecting Xplor vs. PDBv2.3 vs. PDBv3.0 atom" << endl; 
   cerr  << "                       names in reduce.cpp and flipmemo.cpp. Updated the Version and shortVersion" << endl; 
   cerr  << "7/31/07 - rmi -        Bug fix: changed the logic in selecting Xplor vs. PDBv2.3 vs. PDBv3.0 atom" << endl; 
   cerr  << "                       names in reduce.cpp and flipmemo.cpp. Updated the Version and shortVersion" << endl; 
   cerr  << "                       strings and changed the year to the four digit year in the versionString." << endl; 
   cerr  << "8/01/07 - rmi          Added the -Version and -Changes flags. The -Changes flag uses a new function" << endl; 
   cerr  << "                       reduceChanges and follows the format of reduceHelp" << endl; 
   cerr  << "                       Also commented out non-RNA non-DNA 'nucleic acids' from StdResH.cpp to fix" << endl;
   cerr  << "                       doubling of backbone hydrogens." << endl;
   cerr  << "8/14/07 - rmi          Fixed bug in StdResH.cpp to add hydrogens on MSE. The change was to add SE" << endl;
   cerr  << "                       as well as SED as the name for the selenium atom" <<endl; 
   cerr  << "8/18/07 - rwgk         Patched Elementinfo.cpp for compiler problems: (a)Visual C++ warning and (b)Tru64 error" << endl;
   cerr  << "svn rev 67, 68         (a)threw runtime error on fixAtomName() (b)added 'using std::sprintf'" << endl;
   cerr  << "8/29/07 - rmi          Modified the reduce het dict so that hydrogens are not built on carboxylates" << endl; 
//   cerr  << "9/25/07 - rmi          Added a flag ADDNHATGAP which allows a single hydrogen to be built at the N-termini of chain breaks" << endl; 
   cerr  << "                       added break-amide to StdResH to treat these amides as a special case" << endl; 
   cerr  << "10/3/07 - rmi          Added support for Hybrid36 atom and residue numbers" << endl; 
   cerr  << "11/1/07 - rmi          Added support for two character chainIds" << endl; 
   cerr  << "11/7/07 - rmi          Reverted changes to pdb_sscanf.cpp and changed write format for chains to %-2s" <<endl; 
   cerr  << "11/14/07- rmi          Several BUG fixes:  in PDBrec.h getAtomDescr() wants seqNum not serialNum" <<endl;
   cerr  << "                         in AtomPositions.cpp change format to %-2.2s for two character chains" <<endl;
   cerr  << "                         fixed format strings in read_format.i and write_format.i"  <<endl;
   cerr  << "                         explicitly added Hy36seqNum to pdb++.h" << endl; 
   cerr  << "02/14/08 - vbc & rmi   Bob's fixes for dealing with windows line ending files (no fix for mac files though)." <<endl; 
   cerr  << "             & jjh       I (vbc) attempted to fix some of the warnings for new hydrogen names so molprobity isn't" <<endl; 
   cerr  << "                         quite as swamped." <<endl; 
   cerr  << "02/20/08 - vbc & rmi   Fixed double H bug from Bob's previous correction" << endl;
   cerr  << "02/28/08 - jmw & jjh   Fixed altID bug for H(alpha) when connecting atoms only non-blank altID and only one total conformation" <<endl;
   cerr  << "04/11/08 - jjh         Added -STRING flag to allow scripts in Perl/Python to pass a string to reduce for processing.  Output still directed to standard out." << endl;
   cerr  << "04/28/08 - jjh          fixed 4 character Deuterium recognition w/ PDB 3.0 names" << endl;
   cerr  << "08/21/08 - jjh          added -CHARGEs flag to control charge state output - off by default" << endl;
   cerr  << endl;
   exit(1);
}


// output a list of PDB records
std::ostream& outputRecords(std::ostream& os, const std::list<PDBrec*>& l) {
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
std::istream& inputRecords(std::istream& is, std::list<PDBrec*>& records) {
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
		
		if (active && !drop)
		{
			records.push_back(rec);
		}
		else
		{
			delete rec; rec = 0;
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
      if (strcmp((*(ar.begin()))->chain(), (*(br.begin()))->chain()) != 0) { return FALSE; }

      // cterm C of 'a' must be close to nterm N of 'b'
      // we only look at the first set of conformations

      double gap = distance2((*(ar.begin()))->loc(), (*(br.begin()))->loc());
      if ((1.4-GapWidth) < gap && gap < (1.4+GapWidth)) { return TRUE; } 

      // rmi 070924 add warnings for chain breaks
      else if (gap > (1.4+GapWidth)) { 
         cerr << "*WARNING*: Residues " << (*(ar.begin()))->resname() <<  " " << (*(ar.begin()))->resno() << (*(ar.begin()))->insCode()
              << " and " << (*(br.begin()))->resname() << " " << (*(br.begin()))->resno() << (*(br.begin()))->insCode() << " in chain " 
              << (*(ar.begin()))->chain() << " appear unbonded " << endl << "           "
              << " and will be treated as a chain break" << endl; 
         }
      else if (gap < (1.4-GapWidth)) { 
         cerr << "*WARNING*: Residues " << (*(ar.begin()))->resname() << " " << (*(ar.begin()))->resno() << (*(ar.begin()))->insCode()
              << " and " << (*(br.begin()))->resname() << " " << (*(br.begin()))->resno() << (*(br.begin()))->insCode() << " in chain "
              << (*(ar.begin()))->chain() << " are too close together " << endl << "           "
              << " and will be treated as a chain break" << endl;
      }
   }
   return FALSE;
}

void analyzeRes(CTab& hetdatabase, ResBlk* pb, ResBlk* cb, ResBlk* nb,
				AtomPositions& xyz, std::list<PDBrec*>& waters, 
				std::vector<std::string>& fixNotes, std::list<PDBrec*>& rlst) {
	if (! (cb && cb->valid(rlst))) { return; } // double check

	// add in the previous and next records

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

	// apl 2007/03/07 -- Bug fix: do not add hydrogens to "N" on non-amino acids
	// prev: 	if (ProcessConnHydOnHets || StdResH::ResXtraInfo().match(resname, "AminoAcid")) {
	//
	if ( StdResH::ResXtraInfo().match(resname, "AminoAcid")) {
		StdResH *hn = NULL;
		if (ctype == CONNECTED_RES) {
			hn = StdResH::HydPlanTbl().get("amide");
		}
		else if (ctype == NTERM_RES) {
			bool ispro = resname.find("PRO") != std::string::npos;
			hn = StdResH::HydPlanTbl().get(ispro?"nt-pro":"nt-amide");
		}
		else if (ctype == FRAGMENT_RES) { // don't create NH
                        if (NeutralTermini) {
                           hn = StdResH::HydPlanTbl().get("break-amide"); // rmi 070925
                        }
			if (StandardizeRHBondLengths) {
				// case A - really is a N terminal
				bool ispro = resname.find("PRO") != std::string::npos;
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
		if (hn && hn->exclude().find(resname.c_str()) == std::string::npos) {
			std::list<atomPlacementPlan*> temp = hn->plans();
			app.splice(app.end(), temp);
		}
		if (StdResH::ResXtraInfo().match(resname, "AminoAcid")) {
			StdResH *ha = StdResH::HydPlanTbl().get("alpha");
			if (ha && ha->exclude().find(resname.c_str()) == std::string::npos) {
				std::list<atomPlacementPlan*> temp = ha->plans();
				app.splice(app.end(), temp);
			}
		}
	}

	if (StdResH::ResXtraInfo().match(resname, "NucleicAcid")) {
		StdResH *rpbb = StdResH::HydPlanTbl().get("ribose phosphate backbone");
		if (rpbb && rpbb->exclude().find(resname.c_str()) == std::string::npos) {
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
			std::list<atomPlacementPlan*> temp = ct->genHplans(resname.c_str());
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
			if ( (pp.hasFeature(NOTXPLORNAME) &&  (UseXplorNames || UseOldNames) )
				||   (( pp.hasFeature(   XPLORNAME) && ! pp.hasFeature(USEOLDNAMES)) && ! UseXplorNames) ) { 
				return; // keep our naming conventions straight
			}
			if ( (pp.hasFeature(USENEWNAMES) &&   (UseOldNames || UseXplorNames) )
                                ||   (( pp.hasFeature(USEOLDNAMES) && ! pp.hasFeature(   XPLORNAME)) && ! UseOldNames) 
				||   (( pp.hasFeature(USEOLDNAMES) &&   pp.hasFeature(   XPLORNAME)) 
				&&   (! UseOldNames && ! UseXplorNames)) ) {
                                return; // keep our naming conventions straight
                        }
                        if ( (pp.hasFeature(BACKBONEMODEL) &&   ! BackBoneModel) 
                                ||   (pp.hasFeature(   NOTBBMODEL) && ! BackBoneModel) ) {
                                return; // keep our naming conventions straight
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
					maxalt = std::max(maxalt, nconf[i]);
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
					const PDBrec* cnr = rvv[i][std::min(j, nconf[i]-1)];
					loc[i] = cnr->loc(); //apl 7/3/06 -- FIXING PUSH_BACK BUG

					if (j <= nconf[i]-1) {
						char abc = cnr->alt();   // get the alt loc char to use
						if (abc != ' ' && altId == ' ') { // (use only the first...)
							altId = abc;
							occ = cnr->occupancy();
						}
					}
				}
				if (considerNonAlt) { altId = (*(firstAtoms.begin()))->alt(); occ = (*(firstAtoms.begin()))->occupancy(); }

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

					//newHatom->atomno(0);  substituted next call to re-assign atom numbers
                                        newHatom->Hy36Num(0); 
					newHatom->alt(altId);

					newHatom->occupancy(occ);

					newHatom->elemLabel(" H");
					newHatom->chargeLabel("  "); // blank since we don't know
                    
                    //if charges desired, determine and output
                    if(ShowCharges){
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
                    }
					
					if (visableAltConf(*newHatom, DoOnlyAltA)){ // add hyd. only if visible
						// look for cyclization or other modifications to rot. groups

						if ((numConnAtoms > 2) &&
							! okToPlaceHydHere(*newHatom, pp,
							*(rvv[0][std::min(j, nconf[0]-1)]),
							*(rvv[2][std::min(j, nconf[2]-1)]),
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
								*(rvv[0][std::min(j, nconf[0]-1)]),
								*(rvv[1][std::min(j, nconf[1]-1)]),
								*(rvv[2][std::min(j, nconf[2]-1)]),
								TRUE, DemandRotNH3,
								((DemandRotAllMethyls && pp.hasFeature(ROTATEONDEMAND))
								|| (OKProcessMetMe && pp.hasFeature(ROTATEFLAG))) );
						}
						if (DemandFlipAllHNQs && (!resAlts.empty())) {
							xyz.insertFlip(newHatom, resAlts);
						}
					}
					else {
					  //if not added, delete
					  delete newHatom; newHatom = 0; 
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

	if (hptr && hptr->exclude().find(resname) == std::string::npos) {
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
