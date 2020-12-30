// Name: reduce.cpp
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
// Copyright (C) 2020 ReliaSolve
// **************************************************************
//
//  reduceChanges now contains the CHANGELOG or history info
//

#if defined(_MSC_VER)
#pragma warning(disable:4786)
#pragma warning(disable:4305)
#pragma warning(disable:4800)
#endif

const char *versionString =
     "reduce: version 3.10 12/18/2020, Copyright 1997-2016, J. Michael Word; 2020 ReliaSolve";

const char *shortVersion    = "reduce.3.9.201218";
const char *referenceString =
                       "Word, et. al. (1999) J. Mol. Biol. 285, 1735-1747.";
const char *electronicReference = "http://kinemage.biochem.duke.edu";

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
#include "reduce.h"

#define ABANDONED_RC 1

bool Verbose = TRUE;    // do we write processing notes to stdout?
bool KeepConnections          = TRUE;
bool StandardizeRHBondLengths = TRUE;
bool ProcessConnHydOnHets     = TRUE;
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
bool DoOnlyAltA               = FALSE; //jjh changed default 111118
bool OKProcessMetMe           = FALSE; //cjw changed default 160602
bool OKtoAdjust               = TRUE;
bool ShowCliqueTicks          = TRUE;
bool ShowOrientScore          = FALSE;
bool StringInput              = FALSE;
bool ShowCharges              = FALSE;
bool UseNuclearDistances      = FALSE; //jjh 130111
bool UseSEGIDasChain          = FALSE; //jjh 130503
bool ProcessedFirstModel      = FALSE; //jjh 130718
bool RenameFlip               = FALSE; // SJ - 09/25/2015 - flag to specify if the final PDB file will have the coordinates according to
                                       //the rename atoms flip (is the flag is TRUE) or the new rot hinge dock flip (if the flag is FALSE, this is default). Can be set to true by the -renameflip flag in the commandline.
bool GenerateFinalFlip        = FALSE; // SJ - 09/04/2015 to keep track of when scoring and decision of flips finishes and when the final PDB
                                       //coordinates are being generated. This is set to true after all the calculations are done, unless the RenameFlip flag is TRUE

int MaxAromRingDih    = 10;   // max dihedral angle in planarity check for aromatic rings  120724 - Aram

int MinNTermResNo     = 1;   // how high can a resno be for n-term?
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
bool StopBeforeOptimizing = FALSE;  // If true, handle hydrogen drops/adds and then quit
bool AddWaterHydrogens = TRUE;	  // If true, add phantom hydrogens to waters
bool AddOtherHydrogens = TRUE;	  // If true, add hydrogens to non-water atoms
bool RemoveATOMHydrogens = FALSE; // If true, remove hydrogens from ATOM records
bool RemoveOtherHydrogens = FALSE;// If true, remove hydrogens from non-ATOM records

enum ConnType {NTERM_RES, CONNECTED_RES, FRAGMENT_RES};

SummaryStats Tally;

std::ostream& outputRecords(std::ostream& os, const std::list< std::shared_ptr<PDBrec> >& records, int model); // SJ 08/04/2015 added last argument to keep track of how many models have been printed.
void invalidateRecords(std::list< std::shared_ptr<PDBrec> >& rlst);
void renumberAndReconnect(std::list< std::shared_ptr<PDBrec> >& rlst);
void renumberAtoms(std::list< std::shared_ptr<PDBrec> >& rlst);
void analyzeRes(CTab& db, ResBlk* pb, ResBlk* cb, ResBlk* nb,
				AtomPositions& xyz, std::list< std::shared_ptr<PDBrec> >& waters,
				std::vector<std::string>& fixNotes, std::list< std::shared_ptr<PDBrec> >& rlst);
bool isConnected(ResBlk* a, ResBlk* b);
void genHydrogens(const atomPlacementPlan& pp, ResBlk& theRes, bool o2prime,
				  AtomPositions& xyz, std::list<char>& resAlts,
		  std::vector<std::string>& fixNotes, std::list< std::shared_ptr<PDBrec> >& rlst);
void noteWaterInfo(ResBlk& r, std::list< std::shared_ptr<PDBrec> >& waters);
void findAndStandardize(const char* name, const char* resname, ResBlk& theRes);
void stdBondLen(float dist, PDBrec& ourHydrogen, std::list< std::shared_ptr<PDBrec> >& firstAtoms,
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
   std::list< std::shared_ptr<PDBrec> >& nearr, const char * msg);

int optimize(AtomPositions& xyz, std::vector<std::string>& adjNotes) {
    int ret = 0;

    GenerateFinalFlip = FALSE; // SJ 09/25/2015 - this has to be reset to FALSE, because for a new model the scoring and decision for flips has to be done using renaming the atoms and not the three step flip.

		if (Verbose) {
			if (!OFile.empty()) {
				cerr << "Using orientation info in \"" << OFile << "\"." << endl;
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
				cerr << "Orientation penalty scale = " << PenaltyMagnitude
					<< " (" << int(PenaltyMagnitude * 100.0 + 0.5) << "%)" << endl;
			}
			else {
				cerr << "No Orientation penalty (-penalty0.0)" << endl;
			}
			cerr << "Eliminate contacts within " << NBondCutoff << " bonds." << endl;
			cerr << "Ignore atoms with |occupancy| <= "
				<< OccupancyCutoff << " during adjustments." << endl;
			cerr << "Waters ignored if B-Factor >= " << WaterBcutoff
				<< " or |occupancy| < " << WaterOCCcutoff << endl;
#ifdef AROMATICS_ACCEPT_HBONDS
			cerr << "Aromatic rings in amino acids accept hydrogen bonds." << endl;
#endif
			if (DemandFlipAllHNQs) {
				cerr << "Flipping Asn, Gln and His groups." << endl;
				cerr << "For each flip state, bumps where gap is more than "
					<< BadBumpGapCut << "A are indicated with \'!\'." << endl;

			}
			if (DemandRotNH3) {
				cerr << "Rotating NH3 Hydrogens." << endl;
			}
			if (DemandRotAllMethyls) {
				cerr << "Rotating ALL methyls." << endl;
			}
			if (OKProcessMetMe) {
				if (!DemandRotAllMethyls) {
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
		if (nscnt >= 0) { Tally._num_adj += nscnt; }
	    else { // failed to initialize or too many permutations, make note
	      ret = ABANDONED_RC;
	    }
	 }

// record clique and singleton adjustments
     
     if(!RenameFlip)
         GenerateFinalFlip=TRUE; // SJ 09/25/2015 - If the RenameFlip flag is False, that means that the coordinates in the PDB output should be with the rot hinge dock flip. Therefre this flag has to be set to be TRUE, so that the FlipMemo::setOrientation function does the three step flip now.

     // SJ - 09/25/2015 - have changed the formatClique and formatSingles function to call FlipMemo::setOrientations to do the or hinge dock flip if the GenerateFinalFlip flag is TRUE
	 for (int jj = 0; jj < clst.numCliques(); jj++) {
	    clst.formatClique(adjNotes, jj, xyz); // SJ - added the last argument, as this is needed to do the flip
	 }
	 clst.formatSingles(adjNotes, xyz); // SJ - added the last argument, as this is needed to do the flip
      }

   if (Tally._H_added || Tally._H_removed) {
      ;//renumberAndReconnect(records);
   }

  // std::for_each(records.begin(), records.end(), DeleteObject());
  return ret;
}

// SJ 08/03/2015 for printing all records together
void outputRecords_all(const std::vector <std::list< std::shared_ptr<PDBrec> > >& l, std::ostream& os) {
    
    int model=0; // keeping track of how many models are printed
    for (std::vector<std::list< std::shared_ptr<PDBrec> > >::const_iterator ptr = l.begin(); ptr != l.end(); ++ptr) {
        model++;
        outputRecords(os,(const std::list< std::shared_ptr<PDBrec> > &)*ptr, model); // print
    }
    outputRecords(os,(const std::list< std::shared_ptr<PDBrec> > &)l.back(), 0); // SJ - so that all records after the last ENDMDL are printed. model will be equal to 0, which means all models are printed and only the left over information needs to be printed. model = 0 is a little counterintuitive, but that is the only number I can think of that will work.
    
    return;
}

std::string outputRecords_all_string(const std::vector <std::list< std::shared_ptr<PDBrec> > >& all_records) {
	std::ostringstream ss;
	outputRecords_all(all_records, ss);
	return ss.str();
}

// output a list of PDB records
std::ostream& outputRecords(std::ostream& os, const std::list< std::shared_ptr<PDBrec> >& l, int model) {
    
    bool flag; // SJ 08/04/2015 to keep track of if this record has to be printed or not.
    
    if (model == 1) { // this model does not necessarily correspond to the model number specified in the PDB file. This is just to keep internal records of how many models have been printed. See the for loop in the function above.
        flag=true; // everything has to be printed if this is the first model
        
    //SJ: 08/03/2015 copied over from processPDBFile - to be printed only once
	os << "USER  MOD " << shortVersion;
	if(RemoveATOMHydrogens || RemoveOtherHydrogens)
        {
	  os <<" removed " << Tally._H_removed
            << " hydrogens (" << Tally._H_HET_removed << " hets)";
        }
        if(AddWaterHydrogens || AddOtherHydrogens)
        {
          os  <<" H: found="<<Tally._H_found
            <<", std="<< Tally._H_standardized
            << ", add=" << Tally._H_added
	    << ", rem=" << Tally._H_removed
            << ", adj=" << Tally._num_adj;
        }
	os << endl;
	if (!StopBeforeOptimizing)
	{
	  if (Tally._num_renamed > 0) {
	    os << "USER MOD renamed " << Tally._num_renamed
	      << " ambiguous 'A' atoms (marked in segID field)" << endl;
	  }

	  if (UseSEGIDtoChainMap) {
	    PDBrec::DumpSEGIDtoChainMap(os, "USER  MOD mapped ");
	  }
	}
    }
    else
        flag=false; //for all other models, print only what is within MODEL and ENDMDL
    
    for (std::list< std::shared_ptr<PDBrec> >::const_iterator ptr = l.begin(); ptr != l.end(); ++ptr) {
        
        if((*ptr)->type() == PDB::MODEL){
            if(model != 0) // to make sure the last model is not printed again when the function is called to print the leftover stuff.
               flag=true; // start printing (if not model == 1, in which case flag is already true)
        }
        
		if ((*ptr)->valid() && flag) //print only if flag=true
			os << (const PDBrec&)(**ptr) << endl;
        
        if ((*ptr)->type() == PDB::ENDMDL) {
            if(model != 0)
               flag=false;//stop printing anything else
            else
                flag=true; // to make sure the left over stuff is printed when model = 0
        }
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

// check input records for SEGIDs only, use instead of chainID
bool checkSEGIDs(std::list< std::shared_ptr<PDBrec> >& rlst) {
  bool ret = FALSE;
  int full_chain_ctr = 0;
  int full_segid_ctr = 0;
  typedef std::list< std::shared_ptr<PDBrec> >::iterator pdb_iter;
  for (pdb_iter it = rlst.begin(); it != rlst.end(); ) {
	 std::shared_ptr<PDBrec> r = *it;
    //currently only checks atom records - is this enough?
    if (r->type() == PDB::ATOM) {
      if (strcmp(r->chain(), "") != 0) {
        full_chain_ctr++;
      }
      if (strcmp(r->segidLabel(), "") != 0) {
        full_segid_ctr++;
      }
    }
    it++;
  }
  //cerr << "full_chain_ctr = " << full_chain_ctr << endl;
  //cerr << "full_segid_ctr = " << full_segid_ctr << endl;
  if ( (full_chain_ctr == 0) && (full_segid_ctr > 0) ) {
    //cerr << "Using SEGID as chain" << endl;
    ret = TRUE;
  }
  return ret;
}

// input a list of PDB records
std::list< std::shared_ptr<PDBrec> > inputRecords(std::istream& is, int &ModelToProcess, int &ModelNext, int &ModelActive) {
  std::list< std::shared_ptr<PDBrec> > records;
  PDB inputbuffer;
  bool active = TRUE;
  bool modelactive = FALSE;  //041113

  while ((is >> inputbuffer).gcount() > 0) {
	std::shared_ptr<PDBrec> rec = std::make_shared<PDBrec>(inputbuffer);
    bool drop = FALSE;

    switch (rec->type()) {
      case PDB::MODEL:
        if (!ProcessedFirstModel) {
          ModelToProcess = rec->modelNum();
          ProcessedFirstModel = TRUE;
        }
        if (ModelToProcess == rec->modelNum()) {
          active = TRUE;
          modelactive = TRUE; //041113
        }
        else {
          active = FALSE;
          if(modelactive) {//041113
            modelactive = FALSE;
	    /*next after one just done*/
	    ModelNext = rec->modelNum();
          }
          if (Verbose) {
            //cerr << "NOTE: skipping model " << rec->modelNum() << endl;
            cerr << "Processing Model " << ModelToProcess << ", for now skipping model " << rec->modelNum() << endl; //041113
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
  return records;
}

std::vector< std::list< std::shared_ptr<PDBrec> > > inputModels(std::string s)
{
  int ModelToProcess = 1;
  int ModelNext = 0;
  int ModelActive = 0;
  std::vector< std::list< std::shared_ptr<PDBrec> > > models;
  while (ModelToProcess) {
    std::stringstream ss(s);
    models.push_back(inputRecords(ss, ModelToProcess, ModelNext, ModelActive));
    if (ModelNext > 0) {
      ModelToProcess = ModelNext;
      ModelNext = 0; /*perhaps to be rediscovered in PDB file*/
    } else {
      /*ModelNext==0, time to quit*/
      ModelToProcess = 0;
    }
  }

  return models;
}

void renumberAndReconnect(std::list< std::shared_ptr<PDBrec> >& rlst) {
	std::list< std::shared_ptr<PDBrec> >::iterator it = rlst.begin();

	if (KeepConnections && Tally._conect > 0) {
		std::list< std::shared_ptr<PDBrec> > conn;
		std::map<long, std::shared_ptr<PDBrec> > atomsBySeqNum;

		// first we organize atom records by the original sequence num
		// and put all the connect records in a list

		for (; it != rlst.end(); ++it) {
			std::shared_ptr<PDBrec> rin = *it;
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

		for (std::list< std::shared_ptr<PDBrec> >::iterator ic = conn.begin(); ic != conn.end(); ++ic) {
			std::shared_ptr<PDBrec> rup = *ic;
			int anum[11];
			rup->getConect(anum);

			for (int k=0; k < 11; k++) {
				if (anum[k] > 0) {
					std::map<long, std::shared_ptr<PDBrec> >::const_iterator iter = atomsBySeqNum.find(anum[k]);
					if (iter != atomsBySeqNum.end())
						anum[k] = iter->second->atomno();
					else
						anum[k] = -999;
					//	       std::shared_ptr<PDBrec> x = atomsBySeqNum.get(anum[k]);
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

void invalidateRecords(std::list< std::shared_ptr<PDBrec> >& rlst) {
  typedef std::list< std::shared_ptr<PDBrec> >::iterator pdb_iter;
  for (pdb_iter it = rlst.begin(); it != rlst.end(); ) {
	std::shared_ptr<PDBrec> r = *it;
    if (r->type() == PDB::ATOM) {
	  if (r->isHydrogen()) {
        Tally._H_removed++;
      }
    }
    else if (r->type() == PDB::HETATM) {
      if (r->isHydrogen()) {
        Tally._H_removed++;
        Tally._H_HET_removed++;
      }
    }
    r->invalidateRecord();
    it++;
  }
}

void dropHydrogens(std::list< std::shared_ptr<PDBrec> >& rlst, bool RemoveATOMHydrogens, bool RemoveOtherHydrogens) {
  typedef std::list< std::shared_ptr<PDBrec> >::iterator pdb_iter;
	for (pdb_iter it = rlst.begin(); it != rlst.end(); ) {
		std::shared_ptr<PDBrec> r = *it;
    bool need_increment = true;
		if (RemoveATOMHydrogens && (r->type() == PDB::ATOM)) {
			if (r->isHydrogen()) {
				Tally._H_removed++;
				it = rlst.erase(it);
				need_increment = false;
			}
			Tally._num_atoms++;
		}
		else if (RemoveOtherHydrogens) {
			if (r->type() == PDB::HETATM) {
				if (r->isHydrogen()) {
					Tally._H_removed++;
					Tally._H_HET_removed++;
					it = rlst.erase(it);
					need_increment = false;
				}
				Tally._num_atoms++;
			}
			else if (r->type() == PDB::SIGATM
				|| r->type() == PDB::ANISOU
				|| r->type() == PDB::SIGUIJ) { // supplemental records
				if (r->isHydrogen()) {
					it = rlst.erase(it);
					need_increment = false;
				}
			}
			else if (r->type() == PDB::CONECT) {
				Tally._conect++;
			}
		}
		if (need_increment) it++;
	}
  if (Verbose) {
    cerr << "Trimming: removed " << Tally._H_removed << " hydrogens ("
      << Tally._H_HET_removed << " hets)" << endl;
  }
}

// count record types and group records by positon
// also, update iterator to point to the start of the atom records

void scanAndGroupRecords(std::list< std::shared_ptr<PDBrec> >& rlst, AtomPositions& xyz,
						 std::list< std::shared_ptr<PDBrec> >::iterator& startAtoms) {
	bool foundStart = FALSE;

	for (std::list< std::shared_ptr<PDBrec> >::iterator it = rlst.begin(); it != rlst.end(); ++it) {
		std::shared_ptr<PDBrec> r = *it;
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

void reduceList(CTab& hetdatabase, std::list< std::shared_ptr<PDBrec> >& rlst,
				AtomPositions& xyz, std::vector<std::string>& fixNotes) {
	std::list< std::shared_ptr<PDBrec> >::iterator it = rlst.begin();
	std::list< std::shared_ptr<PDBrec> >::iterator it_backup;
	ResBlk *pb = NULL, *cb = NULL, *nb = NULL;
	std::list< std::shared_ptr<PDBrec> > waters;

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

	if (AddWaterHydrogens) { xyz.generateWaterPhantomHs(waters); }
}

void renumberAtoms(std::list< std::shared_ptr<PDBrec> >& rlst) {
	int ia = 0;

	for (std::list< std::shared_ptr<PDBrec> >::iterator it = rlst.begin(); it != rlst.end(); ++it) {
		std::shared_ptr<PDBrec> r = *it;
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

   std::list< std::shared_ptr<PDBrec> > ar;
   a->get(" C", ar);
   std::list< std::shared_ptr<PDBrec> > br;
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

// is the hydrogen atom in methyl group stemmed from aromatic ring? - Aram 07/18/12
bool isAromMethyl(ResConn& ct, const atomPlacementPlan& pp, ResBlk& theRes, const char* resname) {

	std::list<std::string> temp = ct.findRingBondedToMethyl(pp.name(), resname);
	int AROMATIC_RING_DIHEDRAL = MaxAromRingDih;

	if (temp.empty()) {
		return FALSE;
	}

	std::list<std::string>::const_iterator it = temp.begin();
	std::vector< std::shared_ptr<PDBrec> > r_vec;
	r_vec.reserve(temp.size());

	// std::cout << std::endl << "resname: " << resname << ", atomname: " << pp.name() << "." << theRes.firstRec().atomname();
	while (it != temp.end()) {
		std::list< std::shared_ptr<PDBrec> > r_list;
		// if atom does not exist, do not test dihedral
		if (!theRes.contains(*it)) {
		  cerr <<"WARNING: No " << *it << " atom! Cannot determine if "
		       << pp.name() << " is part of an aromatic methyl." << endl;
		  return FALSE;
		}
		theRes.get(*it, r_list);
		std::shared_ptr<PDBrec> rec = *r_list.begin(); // Do not consider alt configulation for now
		r_vec.push_back(rec);
		// std::cout << " " << *it << "(" << r_vec[r_vec.size()-1]->loc() << ")";
		++it;
	}
	// std::cout << std::endl;

	for (int i = 0; i < r_vec.size(); i++) { // check dihedral angles starting from each atom in the ring
		int i0=i, i1=(i+1)%r_vec.size(), i2=(i+2)%r_vec.size(), i3=(i+3)%r_vec.size();
		//cerr << r_vec[i0]->atomname() << r_vec[i1]->atomname() << r_vec[i2]->atomname() << r_vec[i3]->atomname() << endl;
		//cerr << r_vec[i0]->loc() << r_vec[i1]->loc() << r_vec[i2]->loc() << r_vec[i3]->loc() << endl;
		double dih = dihedral(r_vec[i0]->loc(), r_vec[i1]->loc(), r_vec[i2]->loc(), r_vec[i3]->loc());
		//need to check for this? - JJH 130326
		/*if (std::isnan(dih)) {
		  cerr << "problem with dihedral" << endl;
		}*/
		while (dih < 0) dih += 360;
		while (dih >= 180) dih -= 180;
		if ( (dih > AROMATIC_RING_DIHEDRAL) && (dih < (180.0 - AROMATIC_RING_DIHEDRAL)) ) {
			// std::cout << "    dihedral test failed: " << r_vec[i]->atomname() << " (" << dih << ")" << std::endl;
			return FALSE;
		} else {
			// std::cout << ", " << dih;
		}
	}

	return TRUE;
}

//SJ - called for each residue
void analyzeRes(CTab& hetdatabase, ResBlk* pb, ResBlk* cb, ResBlk* nb,
				AtomPositions& xyz, std::list< std::shared_ptr<PDBrec> >& waters,
				std::vector<std::string>& fixNotes, std::list< std::shared_ptr<PDBrec> >& rlst) {
	if (! (cb && cb->valid(rlst))) { return; } // double check

	// add in the previous and next records

	if (pb) {
		std::list< std::shared_ptr<PDBrec> > pr_list;
		pb->get(" C", pr_list);
		for (std::list< std::shared_ptr<PDBrec> >::iterator pr = pr_list.begin(); pr != pr_list.end(); ++pr) {
			cb->addPrevRec(*pr);
		}
	}
	if (nb) {
		std::list< std::shared_ptr<PDBrec> > nr_list;
		nb->get(" N", nr_list);
		for (std::list< std::shared_ptr<PDBrec> >::iterator nr = nr_list.begin(); nr != nr_list.end(); ++nr) {
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
	std::list<std::shared_ptr<atomPlacementPlan> > app;

	// apl 2007/03/07 -- Bug fix: do not add hydrogens to "N" on non-amino acids
	// prev: 	if (ProcessConnHydOnHets || StdResH::ResXtraInfo().match(resname, "AminoAcid")) {
	//
	if ( StdResH::ResXtraInfo().match(resname, "AminoAcid")) { // S.J. this if statement gets the plan for H atoms for the N atom and the Calpha atom given in the StdResH.cpp
        StdResH *hn = NULL;
		if (ctype == CONNECTED_RES) {
			hn = StdResH::HydPlanTbl().get("amide");
		}
		else if (ctype == NTERM_RES) {
			bool ispro = resname.find("PRO") != std::string::npos;
			hn = StdResH::HydPlanTbl().get(ispro?"nt-pro":"nt-amide");
			//Trim N-H atom if present in input model for N-terminal amino JJH 130110
			std::list< std::shared_ptr<PDBrec> > cr_list;
		    cb->get(" H", cr_list);
		    if (cr_list.size() > 0) {
		      std::cerr << "Removing redundant N-terminal N-H hydrogen atom from input model." << std::endl;
		      invalidateRecords(cr_list);
		    }
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
			std::list< std::shared_ptr<PDBrec> > nfr_list;
			cb->get(" N", nfr_list);
			for (std::list< std::shared_ptr<PDBrec> >::iterator nfr = nfr_list.begin(); nfr != nfr_list.end(); ++nfr) {
				(*nfr)->elem(* ElementInfo::StdElemTbl().element("Nacc"));
			}
		}
		if (hn && hn->exclude().find(resname.c_str()) == std::string::npos) {
			std::list<std::shared_ptr<atomPlacementPlan> > temp = hn->plans();
			app.splice(app.end(), temp);
		}
		if (StdResH::ResXtraInfo().match(resname, "AminoAcid")) {
			StdResH *ha = StdResH::HydPlanTbl().get("alpha");
			if (ha && ha->exclude().find(resname.c_str()) == std::string::npos) {
				std::list<std::shared_ptr<atomPlacementPlan> > temp = ha->plans();
				app.splice(app.end(), temp);
			}
		}
	}

	if (StdResH::ResXtraInfo().match(resname, "NucleicAcid")) { // S.J. this does the same thing as above, but just for Nucleotc acids backbone atoms
		StdResH *rpbb = StdResH::HydPlanTbl().get("ribose phosphate backbone");
		if (rpbb && rpbb->exclude().find(resname.c_str()) == std::string::npos) {
			std::list<std::shared_ptr<atomPlacementPlan> > temp = rpbb->plans();
			app.splice(app.end(), temp);
		}
	}
	bool o2prime = cb->contains(" O2*") || cb->contains(" O2'"); // distinguish RNA and DNA

	// first look in the standard H table...

    StdResH *srh = StdResH::HydPlanTbl().get(resname);
	if (srh) {
        std::list<std::shared_ptr<atomPlacementPlan> > temp = srh->plans(); // S.J. this has the hydrogen atom plans for this specific residue given in StdResH.cpp

		// for aromatic methyls - Aram 07/23/12
		for(std::list<std::shared_ptr<atomPlacementPlan> >::const_iterator iter = temp.begin(); iter != temp.end(); ++iter) {
            
            std::string conn2atom = (*iter)->conn(1);
			bool Arom = StdResH::ResXtraInfo().atomHasAttrib(resname, conn2atom, AROMATICFLAG);
			if (Arom) (*iter)->addFeature(AROMATICFLAG);
		}

		app.splice(app.end(), temp);
	}

	else if (ProcessConnHydOnHets) { // if not in the std table we look in the het database...

		std::shared_ptr<ResConn> ct = hetdatabase.findTable(resname);
		if (ct) {
			std::list<std::shared_ptr<atomPlacementPlan> > temp = ct->genHplans(resname.c_str());

			// for aromatic methyls - Aram 07/18/12
			for(std::list<std::shared_ptr<atomPlacementPlan> >::const_iterator iter = temp.begin(); iter != temp.end(); ++iter) {
				bool AromMethyl = isAromMethyl(*ct, **iter, *cb, resname.c_str());
				if (AromMethyl) (*iter)->addFeature(AROMATICFLAG);
			}

			app.splice(app.end(), temp);
		}
		else {
			cerr << "*WARNING*: Res \""
				<< resname << "\" not in HETATM Connection Database. Hydrogens not added."
				<< endl;
		}
	}

	// work through each atom placement plan - S.J. each atom placement plan is a potential H atom that needs to be added. The hydrogens that are already present are not looked at right now.
	if (AddOtherHydrogens) {
		for (std::list<std::shared_ptr<atomPlacementPlan> >::const_iterator iter = app.begin(); iter != app.end(); ++iter) {
			genHydrogens(**iter, *cb, o2prime, xyz, resAlts, fixNotes, rlst);
		}
	}
}

//SJ - called for each potential hydrogen to be added to the residue
void genHydrogens(const atomPlacementPlan& pp, ResBlk& theRes, bool o2prime,
				  AtomPositions& xyz, std::list<char>& resAlts,
				  std::vector<std::string>& fixNotes, std::list< std::shared_ptr<PDBrec> >& rlst) {

	std::list< std::shared_ptr<PDBrec> > ourHydrogens;
	theRes.get(pp.name(), ourHydrogens); // S.J.- checks if this specific Hydrogen to be added is already present in the residue or not.
	std::list< std::shared_ptr<PDBrec> > firstAtoms;
	theRes.get(pp.conn(0), firstAtoms);
    
    bool doNotAdjustSC = FALSE;
	if (!firstAtoms.empty()) { // must connect to something!

		if (!ourHydrogens.empty()) {// Hydrogen exists
            
			std::shared_ptr<PDBrec> o;
			for (std::list< std::shared_ptr<PDBrec> >::iterator it = ourHydrogens.begin(); it != ourHydrogens.end(); ++it) {
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
					(DemandRotExisting || DemandRotNH3 || pp.hasFeature(AROMATICFLAG))) ) {

					char hac = o->alt();
					PDBrec r0atom, r1atom, r2atom;        // find connecting atoms
					int connatomcount = findConnAtoms(pp, theRes, hac,
						r0atom, r1atom, r2atom);

					if (connatomcount == 3) {
						// for heme methyls - Aram 05/31/12
						if ((DemandRotExisting && pp.hasFeature(ROTATEONDEMAND)
							&& pp.hasFeature(AROMATICFLAG))) {
							xyz.insertRotAromMethyl(*o, r0atom, r1atom, r2atom);
							//std::cout << " in genHydrogens " << o->resname() << pp.name() << std::endl;
						} else {
							xyz.insertRot(*o, r0atom, r1atom, r2atom,
								(DemandRotExisting || RotExistingOH),
								DemandRotNH3,
								(  (DemandRotExisting && DemandRotAllMethyls
								&& pp.hasFeature(ROTATEONDEMAND))
								|| (DemandRotExisting && OKProcessMetMe
								&& pp.hasFeature(ROTATEFLAG))) );
						}

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
            int i = 0, j = 0, k = 0;

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
			std::vector<char> all_confs;
			nconf.resize(numConnAtoms);

			std::vector< std::vector< std::shared_ptr<PDBrec> > > rvv;
			rvv.reserve(numConnAtoms);
			bool success = TRUE, alt_success = TRUE;

			for (i = 0; i < numConnAtoms; i++) { // get the records and count how many
				std::list< std::shared_ptr<PDBrec> > rs;
				theRes.get(pp.conn(i), rs);
				if (!rs.empty()) {
					nconf[i] = rs.size();
					maxalt = std::max(maxalt, nconf[i]);
					std::vector< std::shared_ptr<PDBrec> > rvv_v;
					rvv_v.reserve(nconf[i]);
					std::list< std::shared_ptr<PDBrec> >::iterator it_rs = rs.begin();
					for(j=0; j < nconf[i]; ++j, ++it_rs) {
						rvv_v.push_back(*it_rs);
						all_confs.push_back(rvv_v[j]->alt());
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

            // only keep unique chains
            /*(for(k=all_confs.size()-1; k > 0; k--){
                if ( toupper(all_confs[k])
				    < toupper(all_confs[k-1]) ) {
					swap2(all_confs[k], all_confs[k-1]);
			    }
            }*/
            sort (all_confs.begin(), all_confs.end());
            for(k=all_confs.size()-1; k > 0; k--) {
                if ( toupper(all_confs[k])
                    == toupper(all_confs[k-1]) ) {
                    all_confs.erase(all_confs.begin()+k);
                }
            }
            if ( (all_confs.size() > 1) && (all_confs[0] == ' ') ) {
                all_confs.erase(all_confs.begin());
            }

			bool considerNonAlt = FALSE;

			if (pp.hasFeature(STRICTALTFLAG) && (numConnAtoms > 3)
				&& (nconf[0] == 1) && (nconf[1] == 1) && (nconf[2] == 1)) {
				maxalt = 1; // this is an H(alpha) so ignore the fourth (CBeta) alternate conf
				considerNonAlt = TRUE;
				all_confs.clear();
				const std::shared_ptr<PDBrec> cnr = rvv[0][0];
				char abc = cnr->alt();
				all_confs.push_back(abc);
			}

			// LIMITATION:
			// the logic to determine alt conf codes does not handle the case were the chain
			// of atoms switches codes

			std::vector<Point3d> loc(numConnAtoms);
			std::vector<int> counter(numConnAtoms);
			for(j=0; success && j < maxalt; j++) { // for each alternate conformation...
				char altId = ' ', foundId = ' ';
				float occ = (*(firstAtoms.begin()))->occupancy();
                alt_success = TRUE;
				if (considerNonAlt) {
				    const  std::shared_ptr<PDBrec> cnr = rvv[3][j];
				    char abc = cnr->alt();
				    altId = abc;
				}
				else {
				    altId = all_confs[j];
				}
				for(i=0; i < numConnAtoms; i++) {
					const std::shared_ptr<PDBrec> cnr = rvv[i][std::min(j, nconf[i]-1)];
					loc[i] = cnr->loc(); //apl 7/3/06 -- FIXING PUSH_BACK BUG
					counter[i] = std::min(j, nconf[i]-1);
					char abc = cnr->alt();
					if (abc == altId) {
					    if ( (altId != ' ') && (foundId == ' ') ) {
					        foundId = altId;
					        occ = cnr->occupancy();
					    }
					}
					else {
			            alt_success = FALSE;
			            for(k=maxalt-1; k>=0; k--) { //comprehensive search
			                const std::shared_ptr<PDBrec> cnr_v = rvv[i][std::min(k, nconf[i]-1)];
			                char abc = cnr_v->alt();
			                if (abc == altId || abc == ' ') {
				                loc[i] = cnr_v->loc();
				                if ( (abc != ' ') && (foundId == ' ') ) {
				                    foundId = altId;
				                    occ = cnr_v->occupancy();
				                }
				                counter[i] = std::min(k, nconf[i]-1);
				                alt_success = TRUE;
				                break;
			                }
			            }
			            if (!alt_success) {
			                if ( (i>0) && (nconf[i]==1) ) {
			                    const std::shared_ptr<PDBrec> cnr_v = rvv[i][0];
			                    loc[i] = cnr_v->loc();
			                    counter[i] = 0;
				                alt_success = TRUE;
			                }
			                else {
			                    //cerr << "FAILED " << pp.name().substr(0,4).c_str() << " " << altId << endl;
			                    break;
			                }
			            }
					}
				}
				if (considerNonAlt) { altId = (*(firstAtoms.begin()))->alt(); occ = (*(firstAtoms.begin()))->occupancy(); }

				Point3d newHpos;
				if (success && alt_success && (success = pp.placeH(loc, newHpos))) {
					std::shared_ptr<PDBrec> newHatom = std::make_shared<PDBrec>();
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
							*(rvv[0][counter[0]]),
							*(rvv[2][counter[2]]),
							xyz, doNotAdjustSC, fixNotes)) {
							continue; // don't add hyd.
						}

						theRes.insertNewRec(rlst, newHatom);

						noteAddedInTally(*newHatom);

						xyz.put(newHatom); // index in the xyz table

						if (doNotAdjustSC) { continue; } // do not add to the adjustable info

						if ( pp.hasFeature(ROTATEFLAG)
							||  (pp.hasFeature(ROTATEONDEMAND)
							&& (DemandRotAllMethyls || DemandRotNH3 || pp.hasFeature(AROMATICFLAG)) )     ) {
							// for heme methyls - Aram 05/31/12
							if ((pp.hasFeature(ROTATEONDEMAND) && pp.hasFeature(AROMATICFLAG))) {
								//std::cout << " in genHydrogens_noHyd " << newHatom->resname() << pp.name() << std::endl;
								/*xyz.insertRotAromMethyl(*newHatom,
									*(rvv[0][std::min(j, nconf[0]-1)]),
									*(rvv[1][std::min(j, nconf[1]-1)]),
									*(rvv[2][std::min(j, nconf[2]-1)]));*/
								// 130122 - JJH, tracking alternates fix
								xyz.insertRotAromMethyl(*newHatom,
									*(rvv[0][counter[0]]),
									*(rvv[1][counter[1]]),
									*(rvv[2][counter[2]]));
							} else {
								xyz.insertRot(*newHatom,
									/**(rvv[0][std::min(j, nconf[0]-1)]),
									*(rvv[1][std::min(j, nconf[1]-1)]),
									*(rvv[2][std::min(j, nconf[2]-1)]),*/
									// 130122 - JJH, tracking alternates fix
									*(rvv[0][counter[0]]),
									*(rvv[1][counter[1]]),
									*(rvv[2][counter[2]]),
									TRUE, DemandRotNH3,
									((DemandRotAllMethyls && pp.hasFeature(ROTATEONDEMAND))
									|| (OKProcessMetMe && pp.hasFeature(ROTATEFLAG))) );
							}
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
	std::shared_ptr<PDBrec> rec;

	if (pp.num_conn() > 0) {
		std::list< std::shared_ptr<PDBrec> > r0_list;
		theRes.get(pp.conn(0), r0_list);
		for (std::list< std::shared_ptr<PDBrec> >::iterator r0 = r0_list.begin(); r0 != r0_list.end(); ++r0) {
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
		std::list< std::shared_ptr<PDBrec> > r1_list;
		theRes.get(pp.conn(1), r1_list);
		for (std::list< std::shared_ptr<PDBrec> >::iterator r1 = r1_list.begin(); r1 != r1_list.end(); ++r1) {
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
		std::list< std::shared_ptr<PDBrec> > r2_list;
		theRes.get(pp.conn(2), r2_list);
		for (std::list< std::shared_ptr<PDBrec> >::iterator r2 = r2_list.begin(); r2 != r2_list.end(); ++r2) {
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

	std::list< std::shared_ptr<PDBrec> > emptySet;

	const Point3d& theHpos = theHatom.loc();

	if (pp.hasFeature(ROTATEFLAG) || pp.hasFeature(NH3FLAG)) {
		const PDBrec& heavyAtom = a1;
		if (heavyAtom.elem().atno() != 6) { // OH, SH and NH3, but not CH3
			const double halfbondlen = heavyAtom.elem().covRad();
			const double  maxbondlen = halfbondlen
				+ ElementInfo::StdElemTbl().maxCovalentRadius();
			std::list< std::shared_ptr<PDBrec> > nearr_list = xyz.neighbors( heavyAtom.loc(),
				halfbondlen + 0.1,
				maxbondlen  + 0.25);
			std::list< std::shared_ptr<PDBrec> > c2batoms;
			int countOfBonds = 0;
			std::shared_ptr<PDBrec> rec;
			for (std::list< std::shared_ptr<PDBrec> >::const_iterator nearr = nearr_list.begin(); nearr != nearr_list.end(); ++nearr) {
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
								std::shared_ptr<PDBrec> nbs;
								for(std::list< std::shared_ptr<PDBrec> >::const_iterator it = c2batoms.begin(); it != c2batoms.end(); ++it) {
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

		std::list< std::shared_ptr<PDBrec> > nearr_list = xyz.neighbors( nqoxygen.loc(),
			halfbondlen + 0.1, maxbondlen  + 0.25);

		std::shared_ptr<PDBrec> rec;
		for (std::list< std::shared_ptr<PDBrec> >::const_iterator nearr = nearr_list.begin(); nearr != nearr_list.end(); ++nearr) {
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

		std::list< std::shared_ptr<PDBrec> > nearr_list = xyz.neighbors( heavyAtom.loc(),
			halfbondlen + 0.1,
			maxbondlen  + 0.25);
		std::shared_ptr<PDBrec> rec;
		for (std::list< std::shared_ptr<PDBrec> >::const_iterator nearr = nearr_list.begin(); nearr != nearr_list.end(); ++nearr) {
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
		std::list< std::shared_ptr<PDBrec> > nearr_list = xyz.neighbors( theoxygen.loc(),
			pobondlen - 0.25, pobondlen + 0.25);
		std::shared_ptr<PDBrec> rec;
		for (std::list< std::shared_ptr<PDBrec> >::const_iterator nearr = nearr_list.begin(); nearr != nearr_list.end(); ++nearr) {
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
					std::list< std::shared_ptr<PDBrec> >& nearr, const char * msg) {
	std::list< std::shared_ptr<PDBrec> > conAtms;
	std::shared_ptr<PDBrec> rec;
	for(std::list< std::shared_ptr<PDBrec> >::const_iterator it = nearr.begin(); it != nearr.end(); ++it) {
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
		for(std::list< std::shared_ptr<PDBrec> >::const_iterator cat = conAtms.begin(); cat != conAtms.end(); ++cat) {
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
		for(std::list< std::shared_ptr<PDBrec> >::const_iterator it = conAtms.begin(); it != conAtms.end(); ++it) {
			cerr << "-" << (*it)->recName();
		}
		cerr << " " << msg << endl;
	}
}

// For water atoms: ignore any atoms with high b or low occupancy
// and identify waters that do not have explicit H atoms.
// Possible orientations for the new H atoms are identified elsewhere.

void noteWaterInfo(ResBlk& r, std::list< std::shared_ptr<PDBrec> >& waters) {
	std::multimap<std::string, std::shared_ptr<PDBrec> > it_map = r.atomIt();
	std::multimap<std::string, std::shared_ptr<PDBrec> >::const_iterator it = it_map.begin();
	std::string key;
	int i=0, nac=0;
	float occ = 0.0, bf = 0.0;
	const int altbufsz = 10;
	char skipalt[altbufsz], ac;
	bool hwithsamecode;
	std::shared_ptr<PDBrec> a;

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
		std::list<std::shared_ptr<atomPlacementPlan> > app_deque = hptr->plans();
		for (std::list<std::shared_ptr<atomPlacementPlan> >::iterator app = app_deque.begin(); app != app_deque.end(); ++app) {

			const std::shared_ptr<atomPlacementPlan> pp = *app;
			std::list< std::shared_ptr<PDBrec> > firstAtoms;
			theRes.get(pp->conn(0), firstAtoms);

			if (!firstAtoms.empty()) { // must connect to something!

				std::list< std::shared_ptr<PDBrec> > ourHydrogens;
				theRes.get(pp->name(), ourHydrogens);
				for (std::list< std::shared_ptr<PDBrec> >::iterator it = ourHydrogens.begin(); it != ourHydrogens.end(); ++it) {
				  stdBondLen(pp->dist(), **it, firstAtoms, pp->elem());
				}
			}
		}
	}
}

void stdBondLen(float dist, PDBrec& ourHydrogen, std::list< std::shared_ptr<PDBrec> >& firstAtoms,
				const ElementInfo& e) {
	if (ourHydrogen.valid()) {
		std::shared_ptr<PDBrec> temp;
		for(std::list< std::shared_ptr<PDBrec> >::iterator it = firstAtoms.begin(); it != firstAtoms.end(); ++it) {
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
	std::multimap<std::string, std::shared_ptr<PDBrec> > theAtoms_map = theRes.atomIt();
	std::multimap<std::string, std::shared_ptr<PDBrec> >::const_iterator theAtoms = theAtoms_map.begin();
	std::string atomname;
	while (theAtoms != theAtoms_map.end()) {
		atomname = theAtoms->first;
		std::shared_ptr<PDBrec> r;
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
